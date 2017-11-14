# Real-time-QRS-detection
## Purpose
A reliable QRS recognition algorithm would be of critical use in clinical senarios. The code, inspired by Pan & Tompkins, can hopefully realize a real-time analysis of heart rate, Bradycardia, Tachycardia, Premature ventricular contractions and Atrial premature beats.
## Reference
The algorithm is based on the paper: Pan J, Tompkins W J. A real-time QRS detection algorithm[J]. IEEE transactions on biomedical engineering, 1985 (3): 230-236.
## Flow Chart
![Alt text](/imgs/ECG_detect.png)
## Analysis Criteria
![Alt text](/imgs/criteria.png)
## Code Structure
ECG_demo: main function. It will plot the ECG signal in real time.  
rate_cal: calculate real-time heart rate. It has been optimized for time efficiency.  
ECG_diagnosis: abnormality judgement. All the results would be printed in the cammand line.  
## Example
The data has already been uploaded. Say we use 7.txt for evaluation. The results would be as follows:
![Alt text](/imgs/result_screenshot.png)  
After loading all the data, the average heart rate and diagnosis results would be available.
![Alt text](/imgs/result_screenshot2.png)  
All the QRS detection results would be visualized as follows:
![Alt text](/imgs/result_QRSdetect.png)
