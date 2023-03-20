
## kit specific parameter:

# our protocol:
#local soloCBstart=25
#local soloCBlen=10
#local soloUMIstart=17
#local soloUMIlen=8
#local whitelist=$META_DIR'whitelist96'  

# alithea:  
local soloCBstart=1
local soloCBlen=14
local soloUMIstart=15
local soloUMIlen=14
local whitelist=$META_DIR'barcodes_96_V5A_brb'  

## general parameter for STAR alignment:
local outSAMmapqUnique=60
local outSAMunmapped="Within"
local soloStrand="Forward"
local quantMode="GeneCounts"
local outBAMsortingThreadN=8
local soloType="CB_UMI_Simple"
local soloUMIdedup="NoDedup 1MM_All"
local soloCellFilter="None"
local soloBarcodeReadLength=0
local soloFeatures="Gene GeneFull SJ Velocyto"
local outSAMattributes="NH HI nM AS CR UR CB UB GX GN sS sQ sM"
local outFilterMultimapNmax=1
local outSAMtype="BAM SortedByCoordinate"
local outSAMorder="Paired"
local limitBAMsortRAM=60413900847
local readFilesCommand="gunzip -c"

