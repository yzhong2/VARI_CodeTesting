# VARI_CodeTesting  
VARI  code testing for job application.   	
Input data file : coding_challenge_final.vcf      
Project description file : VARI_BRS_position_coding_challenge_instructions_v2.pdf    
Output file : coding_challenge_final_Annot.txt        
R main function file : vcfMain.R    
R file contains all functions : vcfFunc.R   
# Prerequisites 
The program will run in a R enviroment.   
R version : 3.4.0   
Platform: x86_64-pc-linux-gnu (64-bit)   

# Installing    
The program needs several R packages to be installed:    
'readr', 'GenomicAlignments', 'httr', 'jsonlite' and 'tidyr' installed from 
Cran or bioconductor repositories.

# Running the tests
Open "vcfMain.R" file and modify the variable 'fName' into the location of the file coding_challenge_final.vcf in your system.
change the variable 'oFile' into the file name that you want to name for the output of the program in your system.
Change the variable 'odir' into the directory name where you want to put the output file in your system.
Open one R section to run the commands in vcfMain.R using interactive mode after above modification

# Authors
Yi Zhong. Ph.D
