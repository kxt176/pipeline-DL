This repository provides a one-click, PowerShell-driven NGS analysis pipeline for Illumina data on Windows using WSL2. The pipeline automates:
• Reference genome download & indexing
• Quality control (FastQC)
• Read alignment (BWA + SAMtools)
• Post-alignment processing (Picard, GATK)
• Variant calling & compression (GATK → tabix)
• Annotation (SnpEff)
• QC aggregation (MultiQC)
• Summary plots (Python/Matplotlib)

PREREQUISITES

Windows 10 (2004+) or Windows 11

Administrator privileges

Virtualization enabled in BIOS/UEFI

Internet access for downloads

ENABLE VIRTUALIZATION & WSL2
• In BIOS/UEFI, enable virtualization features.
• Open PowerShell as Administrator and run:
dism.exe /online /enable-feature /featurename:VirtualMachinePlatform /all /norestart
dism.exe /online /enable-feature /featurename:Microsoft-Windows-Subsystem-Linux /all /norestart
wsl --set-default-version 2
• Reboot if prompted

INSTALL UBUNTU-22.04 IN WSL2
• In the same elevated PowerShell:
wsl --list --quiet
If “Ubuntu-22.04” is missing, install it:
wsl --install -d Ubuntu-22.04
• Reboot or log out/in, then verify:
wsl --list --quiet

CREATE PIPELINE FOLDER & CLONE
• In an elevated PowerShell:
New-Item -ItemType Directory -Path C:\NGSTools
• Clone the repo into that folder:
git clone https://github.com/your-org/ngs_pipeline.git C:\NGSTools

RUN THE PIPELINE
(Must be an Administrator PowerShell session)
• Execute:
powershell.exe -NoExit -ExecutionPolicy Bypass -STA -File C:\NGSTools\ngs_pipeline_full.ps1
• Follow prompts to select (or skip) a FASTA reference and to choose your FASTQ folder.
• Live progress and logs will display.
• Results are saved under:
C:\NGSTools\results\pipeline_<YYYYMMDD_HHMMSS>

POST-RUN OUTPUTS
In the pipeline_<timestamp> folder you will find:
fastqc_raw/ : FastQC HTML reports
alignment.* : BAM files & metrics
final_variants.vcf.gz + .tbi
final_variants.snpeff.vcf
annotated_variants.tsv
qc_report/ : MultiQC summary
qual_boxplot.png, dp_histogram.png, volcano.png

TROUBLESHOOTING
• “Please run this script as Administrator.” → Relaunch PowerShell as Administrator.
• WSL2 errors → Confirm virtualization is enabled; repeat step 1.
• Missing tools in WSL2 → Ensure Ubuntu-22.04 is installed and default WSL version is 2.

ADDITIONAL NOTES
• Helper Bash scripts are auto-generated with LF line endings to avoid Windows/Unix mismatches.
• Update to the latest PowerShell for best compatibility: https://aka.ms/PSWindows
