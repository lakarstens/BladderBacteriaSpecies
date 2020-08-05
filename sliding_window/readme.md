





# Sliding window analysis

Sliding window analysis (SWA) is the same as a *rolling average* or
*moving average*. For a multi-sequence alignment (MSA) of DNA, a fixed
subset is taken from the start of the MSA, and the amount of information
as measured by the Shannon Entropy is averaged for all columns. Then the
subset is shifted forward by one nucleotide. This continues for the
length of the MSA, in the manner of a sliding window. See the picture
below.

![](../../resources/md_files/sm_sliding_window_entropy.png)

The size of the window will produce a different “trace” of the
information in the MSA. Smaller windows will produce jangly graphs,
while very large windows will produce graphs that are almost flat. The
graphs below show this effect on a mock MSA. There are three variable
regions that are 1, 3 and 10 nucleotides wide, and the average Shannon
entropy for the whole MSA is shown by the green line. A window size of
10 shows all three regions fairly well, while a window size of 25 smear
the first two variable regions together. A window size of 60 hides all
variable regions in the MSA.

![](../../resources/md_files/window_sizes.png)

For a MSA of gene sequences like the 16S rRNA gene, the sizes of the
variable regions are fairly well known, and the average size is 75
nucleotides. We tried 20,70,140 and 170 bp window sizes. Running
`weighted_sized_windows.py` calculated the weights for each species in
the MSA, and then ran a sliding window analysis (SWA) for each window
size. For novel multisequence alignments with variable regions of
unknown locations, you have to try a wider
range.

``` shellsession
(base) C:\Users\Carter\Documents\thesis\thesis_git\sliding_window\src_files>python weighted_sized_window.py redo_weight_for_readme.txt
        Creating weighted dictionary... : [##############################] 99%
        writing weighted entropy values to ../processed_data/weighted_ent_08_03_1426/ktw_16s_type_2019-11-25_34_weight_ent.csv
        writing window size 20 to ../processed_data/specific_windows_weighted_08_03_1426
        writing window size 70 to ../processed_data/specific_windows_weighted_08_03_1426
        writing window size 140 to ../processed_data/specific_windows_weighted_08_03_1426
        writing window size 170 to ../processed_data/specific_windows_weighted_08_03_1426
```

As guessed, a window size of 20 generates a jangly graph, while the
window size of 140 smoothes the variable regions too much.

![](readme_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

For this data set, a window size of 70 worked best at demarcating the 9
variable regions of the 16S rRNA gene sequence, shown in the picture
below.

![](readme_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

There are two things of note that are visible in the graph. The first is
there seems to be a bump in entropy between V4 and V5. The standard V4
primer set described by Caporaso spans the V4 and this little bump. The
second is that the V9 region is almost masked by the steep rise in the
entropy at the end of the MSA. The rise in entropy at the end of the MSA
indicates many insertions and deletions of the sequences compared across
the species represented in this data set.

![](readme_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->
