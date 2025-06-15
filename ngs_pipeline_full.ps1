<#
.SYNOPSIS
  One-click NGS pipeline launcher from Elevated PowerShell using WSL2
.DESCRIPTION
  - Enables WSL2 features and installs Ubuntu-22.04
  - Installs bioinformatics tools inside WSL
  - Prompts for reference FASTA; if none provided, downloads GRCh38 primary assembly
  - Indexes the reference FASTA (BWA, Samtools, Picard)
  - Prompts for FASTQ directory
  - Runs FastQC, BWA alignment, Samtools sort/index, bcftools variant calling, MultiQC, and Python summaries
  - Displays live progress and detailed command logging
#>

# Require Administrator
if (-not ([Security.Principal.WindowsPrincipal] [Security.Principal.WindowsIdentity]::GetCurrent()).IsInRole([Security.Principal.WindowsBuiltinRole]::Administrator)) {
    Write-Error "Please run this script as Administrator."
    exit 1
}

# ——————————————————————————————————————————————————————————————————
#  1) WSL setup & Ubuntu install
# ——————————————————————————————————————————————————————————————————

Write-Host "[WSL Setup] Enabling required features..." -ForegroundColor Cyan
dism.exe /online /enable-feature /featurename:VirtualMachinePlatform /all /norestart | Out-Null
dism.exe /online /enable-feature /featurename:Microsoft-Windows-Subsystem-Linux /all /norestart | Out-Null
wsl --set-default-version 2 | Out-Null

# Ensure Ubuntu-22.04 is installed
$distros = wsl --list --quiet
if ($distros -notcontains 'Ubuntu-22.04') {
    Write-Host "[WSL Setup] Installing Ubuntu-22.04..." -ForegroundColor Cyan
    wsl --install -d Ubuntu-22.04
    Write-Host "Reboot your system and re-run this script." -ForegroundColor Yellow
    exit 0
}
# ——————————————————————————————————————————————————————————————————
#  2) Install core tools in WSL
# ——————————————————————————————————————————————————————————————————

Write-Host "[1/9] Installing tools in WSL..." -ForegroundColor Green
wsl bash -lc '
cat << "EOF" > /tmp/install_tools.sh
#!/usr/bin/env bash
set -euo pipefail
sudo DEBIAN_FRONTEND=noninteractive apt-get update -y
sudo DEBIAN_FRONTEND=noninteractive apt-get install -y wget bwa samtools fastqc sra-toolkit picard multiqc bcftools python3-pip python3 dos2unix
pip3 install --user matplotlib numpy
EOF
chmod +x /tmp/install_tools.sh
bash /tmp/install_tools.sh
rm /tmp/install_tools.sh
'
# ——————————————————————————————————————————————————————————————————
#  3) Reference FASTA selection or download
# ——————————————————————————————————————————————————————————————————

Add-Type -AssemblyName System.Windows.Forms
function Select-File($desc, $filter) {
    $ofd = New-Object System.Windows.Forms.OpenFileDialog
    $ofd.Title  = $desc
    $ofd.Filter = $filter
    if ($ofd.ShowDialog() -eq 'OK') { return $ofd.FileName }
    else { return $null }
}

$refWin    = Select-File 'Select reference FASTA (optional; Cancel to auto-download GRCh38)' 'FASTA Files (*.fa;*.fasta;*.fa.gz)|*.fa;*.fasta;*.fa.gz'
$refWinDir = 'C:\NGSTools\reference'
if (!(Test-Path $refWinDir)) { New-Item -ItemType Directory -Path $refWinDir | Out-Null }

if ($refWin) {
    Write-Host "Reference FASTA selected: $refWin" -ForegroundColor Green
    # Only copy if not already in reference folder
    $srcDir = Split-Path $refWin -Parent
    if ($srcDir -ieq $refWinDir) {
        Write-Host "Source file already in reference folder; skipping copy." -ForegroundColor Cyan
    } else {
        Write-Host "Copying user-provided FASTA to $refWinDir…" -ForegroundColor Green
        Copy-Item -Path $refWin -Destination $refWinDir -Force
    }
    # Decompress if gzipped, preserving correct basename
    $leaf = Split-Path $refWin -Leaf
    if ($leaf -like '*.gz') {
        $faName = [IO.Path]::GetFileNameWithoutExtension($leaf)
        Write-Host "Decompressing $leaf to $faName in WSL..." -ForegroundColor Green
        wsl bash -lc "gunzip -c /mnt/c/NGSTools/reference/$leaf > /mnt/c/NGSTools/reference/$faName"
        $refWsl = "/mnt/c/NGSTools/reference/$faName"
    } else {
        $refWsl = "/mnt/c/NGSTools/reference/$leaf"
    }
} else {
    Write-Host "No reference provided; downloading GRCh38 primary assembly…" -ForegroundColor Yellow

# ——————————————————————————————————————————————————————————————————
#  4) Helper scripts in WSL: index_reference, snpeff_download, snpeff_annotate
# ——————————————————————————————————————————————————————————————————
#download_ref.sh
# ——————————————————————————————————————————————————————————————————
    $downloadHelper = @'
#!/usr/bin/env bash
set -euo pipefail

# Ensure ref folder
mkdir -p /mnt/c/NGSTools/reference
cd /mnt/c/NGSTools/reference

echo "Downloading GRCh38 primary assembly..."
curl -# -o Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz \
    ftp://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

echo
echo "Decompressing..."
gunzip -f Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

exit 0
'@

    # Encode, write and chmod inside WSL
    $b64dl = [Convert]::ToBase64String([System.Text.Encoding]::UTF8.GetBytes($downloadHelper))
    wsl bash -lc "echo '$b64dl' | base64 -d > /mnt/c/NGSTools/download_ref.sh && chmod +x /mnt/c/NGSTools/download_ref.sh"

    # Run the helper
    wsl bash -lc "/mnt/c/NGSTools/download_ref.sh"

    # Set the WSL path for downstream steps
    $refWsl = "/mnt/c/NGSTools/reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
}

# ——————————————————————————————————————————————————————————————————
# index_reference.sh
# ——————————————————————————————————————————————————————————————————

$indexHelper = @'
#!/usr/bin/env bash
set -euo pipefail

REF="$1"
DIR="$(dirname "$REF")"
BASE="$(basename "$REF" .fa)"

# 1) BWA index files
if [[ -e "$REF.amb" && -e "$REF.ann" && -e "$REF.bwt" && -e "$REF.pac" && -e "$REF.sa" ]]; then
  echo "BWA index files exist, skipping"
else
  echo "Building BWA indexes..."
  bwa index "$REF"
fi

# 2) Samtools faidx

if [[ -e "$REF.fai" ]]; then
  echo "Samtools FASTA index exists, skipping"
else
  echo "Building Samtools FASTA index..."
  samtools faidx "$REF"
fi

# 3) GATK sequence dictionary
if [[ -e "$DIR/$BASE.dict" ]]; then
  echo "Sequence dictionary exists, skipping"
else
  echo "Creating sequence dictionary with GATK..."
  gatk CreateSequenceDictionary -R "$REF" -O "$DIR/$BASE.dict"
fi

exit 0
'@

$b64 = [Convert]::ToBase64String([System.Text.Encoding]::UTF8.GetBytes($indexHelper))
wsl bash -lc "echo '$b64' | base64 -d > /mnt/c/NGSTools/index_reference.sh && chmod +x /mnt/c/NGSTools/index_reference.sh"

# ——————————————————————————————————————————————————————————————————
# snpeff_download.sh
# ——————————————————————————————————————————————————————————————————

$snpeffDl = @'
#!/usr/bin/env bash
set -euo pipefail

cd /mnt/c/NGSTools/snpEff

# If already present, skip
if [ -d data/GRCh38.99 ]; then
  echo "GRCh38.99 DB exists, skipping"
  exit 0
fi

echo "Downloading snpEff database GRCh38.99…"
java -Xmx4g -jar snpEff.jar download -v GRCh38.99
'@

$b64 = [Convert]::ToBase64String([System.Text.Encoding]::UTF8.GetBytes($snpeffDl))
wsl bash -lc "echo '$b64' | base64 -d > /mnt/c/NGSTools/snpeff_download.sh; chmod +x /mnt/c/NGSTools/snpeff_download.sh"

# ——————————————————————————————————————————————————————————————————
# snpeff_annotate.sh
# ——————————————————————————————————————————————————————————————————

$annotateHelper = @'
#!/usr/bin/env bash
set -euo pipefail

WORKDIR="$1"
INPUT="$WORKDIR/final_variants.vcf.gz"
OUTPUT="$WORKDIR/final_variants.snpeff.vcf"

# make sure output directory exists
mkdir -p "$(dirname "$OUTPUT")"

echo "Running SnpEff annotation for GRCh38.99…"
export _JAVA_OPTIONS='-Xmx8g'
java -jar /mnt/c/NGSTools/snpEff/snpEff.jar GRCh38.99 "$INPUT" > "$OUTPUT"

echo "✔ SnpEff annotation complete: $OUTPUT"
'@

$b64 = [Convert]::ToBase64String([System.Text.Encoding]::UTF8.GetBytes($annotateHelper))
wsl bash -lc "echo '$b64' | base64 -d > /mnt/c/NGSTools/snpeff_annotate.sh && chmod +x /mnt/c/NGSTools/snpeff_annotate.sh"

# ——————————————————————————————————————————————————————————————————
# generate_summary.sh
# ——————————————————————————————————————————————————————————————————


$helper = @'
#!/usr/bin/env bash
set -euo pipefail

# Usage: generate_summary.sh <output-dir-in-WSL>
outdir="$1"

cat << 'PYTHON_BLOCK' > "$outdir/summarize_variants.py"
#!/usr/bin/env python3
import gzip
import matplotlib.pyplot as plt
import numpy as np

quals, dps = [], []
with gzip.open('final_variants.vcf.gz','rt') as fh:
    for line in fh:
        if line.startswith('#'): continue
        cols = line.split('\t')
        q = float(cols[5])
        info = dict(kv.split('=') for kv in cols[7].split(';') if '=' in kv)
        dp = int(info.get('DP',0))
        quals.append(q)
        dps.append(dp)

# QUAL boxplot
plt.figure()
plt.boxplot(quals)
plt.ylabel('QUAL')
plt.savefig('qual_boxplot.png')

# DP histogram
plt.figure()
plt.hist(dps, bins=30, edgecolor='black')
plt.xlabel('DP')
plt.ylabel('Count')
plt.savefig('dp_histogram.png')

# Volcano‐style plot
plt.figure()
plt.scatter(dps, -np.log10(quals), s=5, alpha=0.4)
plt.xlabel('DP')
plt.ylabel('-log10(QUAL)')
plt.savefig('volcano.png')
PYTHON_BLOCK
'@

# Save to Windows location
Set-Content -Path 'C:\NGSTools\generate_summary.sh' -Value $helper -Encoding ASCII

# Make it executable inside WSL
wsl bash -lc "chmod +x /mnt/c/NGSTools/generate_summary.sh"

# Convert to LF and make executable
wsl bash -lc "dos2unix /mnt/c/NGSTools/generate_summary.sh && chmod +x /mnt/c/NGSTools/generate_summary.sh"

# ——————————————————————————————————————————————————————————————————
#  5) FASTQ selection & prepare output
# ——————————————————————————————————————————————————————————————————

Write-Host "Reference FASTA in WSL: $refWsl" -ForegroundColor Green

# Prompt for FASTQ directory
function Select-Folder($desc) {
    $fbd = New-Object System.Windows.Forms.FolderBrowserDialog
    $fbd.Description = $desc
    if ($fbd.ShowDialog() -eq 'OK') { return $fbd.SelectedPath } else { exit 1 }
}
$fastqWin = Select-Folder 'Select folder with FASTQ files (*.fastq, *.fastq.gz)'
Write-Host "Selected FASTQ folder: $fastqWin" -ForegroundColor Green
$drive, $rest = $fastqWin -split ':',2
$fastqWsl = "/mnt/$($drive.ToLower())$($rest -replace '\\','/')"

# Prepare output folder in WSL
$outBaseWsl = '/mnt/c/NGSTools/results'
wsl bash -lc "mkdir -p $outBaseWsl"
$timestamp = (Get-Date -Format 'yyyyMMdd_HHmmss')
$outDirWsl = "$outBaseWsl/pipeline_$timestamp"
wsl bash -lc "mkdir -p $outDirWsl"

# ——————————————————————————————————————————————————————————————————
#  6) Pipeline steps definition
# ——————————————————————————————————————————————————————————————————

$steps = @{
  Name = 'Download dbSNP'
  Cmd  = @"
cd /mnt/c/NGSTools/reference
if [ -f dbsnp_all.vcf.gz ] && [ -f dbsnp_all.vcf.gz.tbi ]; then
  echo 'dbSNP already exists, skipping'
else
  echo 'Downloading dbSNP (VCF)…'
  curl -# -o dbsnp_all.vcf.gz \
       ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/All_20180418.vcf.gz

  echo
  echo 'Downloading dbSNP index (TBI)…'
  curl -# -o dbsnp_all.vcf.gz.tbi \
       ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/All_20180418.vcf.gz.tbi

  echo
  echo 'dbSNP download complete.'
fi
"@
},
    @{Name='Install GATK'; Cmd="if [ ! -d /mnt/c/NGSTools/reference/gatk-4.6.2.0 ]; then \
    cd /mnt/c/NGSTools/reference && \
    wget -q -O gatk.zip https://github.com/broadinstitute/gatk/releases/download/4.6.2.0/gatk-4.6.2.0.zip && \
    unzip -q gatk.zip && rm gatk.zip; \
  else \
    echo 'GATK already installed, skipping'; \
  fi"},

    @{ Name = 'Install Java (Temurin 21)'
  Cmd  = @"
# 1) If Java 21 is already active, skip the rest
if command -v java >/dev/null 2>&1 && java -version 2>&1 | grep -q 'version \"21'; then
  echo 'Java 21 already installed, skipping'
  exit 0
fi

# 2) Ensure we can add new repos
sudo DEBIAN_FRONTEND=noninteractive apt-get update -y
sudo DEBIAN_FRONTEND=noninteractive apt-get install -y software-properties-common curl

# 3) Add the Adoptium Temurin repo if temurin-21-jre-headless isn’t in cache
if ! apt-cache show temurin-21-jre-headless >/dev/null 2>&1; then
  echo 'Adding Adoptium Temurin 21 repository…'
  curl -fsSL https://packages.adoptium.net/artifactory/api/gpg/key/public \
    | sudo tee /etc/apt/trusted.gpg.d/adoptium.asc >/dev/null
  echo 'deb https://packages.adoptium.net/artifactory/deb jammy main' \
    | sudo tee /etc/apt/sources.list.d/adoptium.list
  sudo DEBIAN_FRONTEND=noninteractive apt-get update -y
fi

# 4) Install Temurin 21, or fall back to OpenJDK 21 if needed
if apt-cache show temurin-21-jre-headless >/dev/null 2>&1; then
  echo 'Installing Temurin 21 JRE…'
  sudo DEBIAN_FRONTEND=noninteractive apt-get install -y temurin-21-jre-headless
else
  echo 'Temurin not found; installing OpenJDK 21…'
  sudo add-apt-repository -y ppa:openjdk-r/ppa
  sudo DEBIAN_FRONTEND=noninteractive apt-get update -y
  sudo DEBIAN_FRONTEND=noninteractive apt-get install -y openjdk-21-jre-headless
fi

# 5) Verify installation
echo 'Installed Java version:'
java -version
"@
},



    @{Name='Index Reference'; Cmd="bash /mnt/c/NGSTools/index_reference.sh '$refWsl'"},
    @{Name='FastQC raw reads';       Cmd="mkdir -p $outDirWsl/fastqc_raw && fastqc -t 4 -o $outDirWsl/fastqc_raw $fastqWsl/*.fastq $fastqWsl/*.fastq.gz"},
    @{Name='Alignment (BWA)';        Cmd="bwa mem -t 4 -R '@RG\tID:sample\tSM:sample\tPL:ILLUMINA' $refWsl $fastqWsl/*.fastq* | samtools view -Sb -o $outDirWsl/alignment.unsorted.bam -"},
    @{Name='Sort BAM';               Cmd="samtools sort -@ 4 -o $outDirWsl/alignment.sorted.bam $outDirWsl/alignment.unsorted.bam"},
    @{Name='Mark Duplicates';        Cmd="java -Xmx4g -jar /mnt/c/NGSTools/reference/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar MarkDuplicates -I $outDirWsl/alignment.sorted.bam -O $outDirWsl/alignment.dedup.bam --METRICS_FILE $outDirWsl/alignment.dupmetrics.txt --CREATE_INDEX true"},
    @{Name='BaseRecalibrator';       Cmd="java -Xmx4g -jar /mnt/c/NGSTools/reference/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar BaseRecalibrator -R $refWsl -I $outDirWsl/alignment.dedup.bam --known-sites /mnt/c/NGSTools/reference/dbsnp_all.vcf.gz -O $outDirWsl/recal.table"},
    @{Name='Apply BQSR';             Cmd="java -Xmx4g -jar /mnt/c/NGSTools/reference/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar ApplyBQSR -R $refWsl -I $outDirWsl/alignment.dedup.bam --bqsr-recal-file $outDirWsl/recal.table -O $outDirWsl/alignment.recal.bam"},
    @{Name='HaplotypeCaller';        Cmd="java -Xmx4g -jar /mnt/c/NGSTools/reference/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar HaplotypeCaller -R $refWsl -I $outDirWsl/alignment.recal.bam -O $outDirWsl/raw_variants.g.vcf.gz -ERC GVCF"},
    @{Name='GenotypeGVCFs';          Cmd="java -Xmx4g -jar /mnt/c/NGSTools/reference/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar GenotypeGVCFs -R $refWsl -V $outDirWsl/raw_variants.g.vcf.gz -O $outDirWsl/final_variants.vcf.gz"},
    @{Name='Index Final VCF';        Cmd="tabix -f -p vcf $outDirWsl/final_variants.vcf.gz"},


@{ Name = 'Install SnpEff'
  Cmd  = @"
# If we've already unpacked snpEff.jar, skip installation
if [ -f /mnt/c/NGSTools/snpEff/snpEff.jar ]; then
  echo 'SnpEff already installed, skipping'
  exit 0
fi

cd /mnt/c/NGSTools

echo 'Installing SnpEff v5.x…'
curl -# -L -o snpEff_latest_core.zip \
     https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip

unzip -q snpEff_latest_core.zip -d /mnt/c/NGSTools
rm snpEff_latest_core.zip

# Verify JAR presence
if [ -f /mnt/c/NGSTools/snpEff/snpEff.jar ]; then
  echo '✔ snpEff.jar found'
else
  echo '❌ ERROR: snpEff.jar not found!' >&2
  exit 1
fi
"@
},

@{ 
  Name = 'Download SnpEff DB'
  Cmd  = "bash /mnt/c/NGSTools/snpeff_download.sh" 
},
  
@{
  Name = 'Annotate (SnpEff)'
  Cmd  = "bash /mnt/c/NGSTools/snpeff_annotate.sh $outDirWsl"
},

  @{ Name = 'VCF to TSV'
     Cmd = 'bcftools query -f ''%CHROM\t%POS\t%REF\t%ALT\t%ANN\n'' ' +
           $outDirWsl + '/final_variants.snpeff.vcf > ' +
           $outDirWsl + '/annotated_variants.tsv' },

  @{ Name = 'MultiQC report'
     Cmd = 'multiqc ' + $outDirWsl + ' -o ' + $outDirWsl + '/qc_report' },

  @{ Name = 'Generate Python summary script'
     Cmd  = "bash /mnt/c/NGSTools/generate_summary.sh `"$outDirWsl`""},

  @{ Name = 'Python summary'
     Cmd  = "cd `"$outDirWsl`" && python3 summarize_variants.py"}

# ——————————————————————————————————————————————————————————————————
#  7) Run with progress bar
# ——————————————————————————————————————————————————————————————————

$total = $steps.Count

$total = $steps.Count
for ($i=0; $i -lt $total; $i++) {
    $step = $steps[$i]
    $percent = [int]((($i)/$total)*100)
    Write-Progress -Activity 'NGS Pipeline' -Status $step.Name -PercentComplete $percent
    Write-Host "`n[Step $($i+1)/$total] $($step.Name)" -ForegroundColor Yellow
    Write-Host "-> $($step.Cmd)" -ForegroundColor Gray
    wsl bash -lc "$($step.Cmd)"
    if ($LASTEXITCODE -ne 0) {
        Write-Error "Step '$($step.Name)' failed with code $LASTEXITCODE"; break
    }
}
Write-Progress -Activity 'NGS Pipeline' -Status 'Finished' -PercentComplete 100 -Completed
Write-Host "Pipeline complete! Results: C:\NGSTools\results\pipeline_$timestamp" -ForegroundColor Green
