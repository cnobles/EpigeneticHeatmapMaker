# EpigeneticHeatmapMaker

#Usage
For list of GTSP and replicates generate heatmap
for a given reference genome and given integration site database:
```
Rscript epi_heatmap_from_db.R sampleName_GTSP.csv -o heatmap --ref_genome hg18  -sites_group intSitesDev237
```

At present we only support "hg18" and "mm8" genomes.

Group should be present in ~/.my.cnf.
