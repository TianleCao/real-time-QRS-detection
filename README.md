# Real-time-QRS-detection
## Purpose
A reliable QRS recognition algorithm would be of critical use in clinical senarios. The code, inspired by Pan & Tompkins, can hopefully realize a real-time analysis of heart rate, Bradycardia, Tachycardia, Premature ventricular contractions and Atrial premature beats.
## Reference
The algorithm is based on the paper: Pan J, Tompkins W J. A real-time QRS detection algorithm[J]. IEEE transactions on biomedical engineering, 1985 (3): 230-236.
## Flow Chart
![Alt text](/path/to/imgs/ECG_detect.png)
## Analysis Criteria
| Category        | Criteria           |
| ------------- |-------------|
|Bradycardia     |$RR_t>1.5s$ or $AR_t>1.2s$
|Tachycardia     |$$AR_t<0.6s$$
|Premature ventricular contractions      |$$RR_{t-1}<0.875AR_{t-2}$$ and $$QRS width>0.12s$$ and $$RR_{t-1}+RR_t=2AR_{t-2}$$
|Atrial premature beats                  |$$RR_{t-1}<0.875AR_{t-2}$$ and $$QRS width normal$$ and $$RR_{t-1}+RR_t<2AR_{t-2}$$
## Code Structure

