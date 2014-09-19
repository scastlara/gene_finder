#############################################################################
####                            GENE_FINDER                          ########
#############################################################################

# *******************************************************************
# These scripts allow you to tag genes in sentences using the 
# HUGO Gene Nomenclature Comittee (HGNC) dictionary
# *******************************************************************

# *******************************************************************
# Every file inside this project is free software, feel free to use it, 
# change it or distribute it.
# ********************************************************************


# *********************************
# Scripts                         #
# *********************************

#**************************************
gene_finder.pl 
  -  This script tags genes using the HUGO dictionary.

        Arguments: 
        - Gene synonyms table (obtained using HUGO_synonyms_downloader.pl)
        - Stopwords file: list of words that will not be analyzed.
          This reduces the number of false positives.
        - Text file(s) you want to analyze.
                  
#**************************************
HUGO_synonyms_downloader.pl 
  - This script downloads a dictionary of gene synonyms
    from HGNC. You can change the source of the synonys, 
    but you'll have to change gene_finder.pl to parse those
    correctly.
                
        
