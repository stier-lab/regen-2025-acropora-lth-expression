Subject: Re: updated LTH summary HTML

Hi Molly,

Thanks for the mock-up. I took another pass at the stage-timing figure with the real data.

I think your layout works well for the basic question: when did each wound-healing stage first show up? The one wrinkle is that the data are scored in whole days, and the wounded subset is 12 fragments per temperature. So the violin version gets a little thin and spiky in places, even though the underlying point is clear.

I made three versions so we can pick the cleanest one:

1. A count-scaled violin plot, closest to your mock-up. It keeps fragments that did not reach a stage in a gray band at the top.

2. A cumulative percent-reached plot. This is probably the clearest for the main summary. It shows that both temperatures get through early healing and tip formation, but 31 °C stalls at new corallite budding. By Day 15, only 4/12 heated fragments reached that final stage.

3. A percent-reached heatmap. This is compact and might work better as a supplement or backup table.

My vote is to use the cumulative plot in the summary, keep the violin version available since it follows your mock-up, and use the heatmap or mean/IQR table as the compact backup.

I also fixed the cohort language in the summary. The main physiology cohort is 48 fragments total. For this stage-timing plot, the wounded subset is 12 fragments per temperature. The microscope-photo cohort is 16 all-wounded fragments total, and I kept that separate for coenosarc healing only. In the current export, 14/14 scored corals had coenosarc cover by Day 1, and 16/16 had it by Day 2.

I updated the same HTML file here:
https://drive.google.com/file/d/18z9zDpKl690ZjBZrhjjIrdtDvllNllPK/view?usp=drivesdk

Best,
Adrian
