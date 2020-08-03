Kidney Transplant Risk Calculator
=======================================

Task
-------------------
DATA3888 Capstone Project Shiny App Only

Mar 2020 - May 2020

Description
-------------------
Kidney failure is the final stage of renal disease and poses a major threat to the body as the excretory system fails to function properly. To combat kidney failure patients can choose two forms of treatment in terms of medical intervention; renal dialysis or organ transplantation. While organ transplantation is greatly preferred, kidney organ allocation has posed itself as a major resource allocation problem. Donors and patients need to be matched effectively and accurately to each other to not only preserve the life of the patient but also, maximise functionality of these limited kidney organs.

With this problem in mind, we developed a tool to aid in the effective and accurate allocation of donor organs to their respective patients. The developed risk calculator will assist practitioners in their decision making and shall even in- form the prescriptions for immunosuppressive drugs. The risk calculator was developed with the intention that it would be used in a clinical setting where shared decision making is implemented. According to Gordan (2013), shared decision making promotes patient centered care. It permits the integration of the nephrologist’s expertise on renal allograft dysfunction with the patient’s values and beliefs concerning future treatment. Within this clinical setting, our hope would be that the risk calculator provides an opportunity of discussion that concerns the nature of treatment prior to, during and post organ transplantation.

Conclusion
-------------------
Through multiple regression and variable selection, a fitted model with solely body measurements was determined for each of the three obesity indicators. Using R2, BMI was identified as the best indicator as it has the highest proportion of variance that can be explained using only body measurements.

Abdomen was the most important body measurement for determining obesity because for all three indicators, it ranked the highest in terms of relative importance in prediction.

By separating the dataset with the binary variable for over-weight individuals, a simplified prediction model with 92% accuracy was built. The simplified model contains three body measurements - chest, abdomen and thigh, and should be relatively simpler to measure.

Challenges & Learning Points
-------------------
1. Analysis cannot be applied to whole population
   - gender bias (only men)

2. Multicollinearity - Variables waist (dropped) and abdomen are highly correlated

3. Using multiple regression to determine which variables to drop in order to attain the final fitted model

4. Using binary logistic regression to interpret relationship between overweight and other significant variables

Files
-------------------
1. DATA2002_Module_4_Report.html
   - Contains main R code and our study findings

2. presentation.html
   - 27 page slides created using xaringan

3. Executive_Report.html
   - A professional 3-page report in the form of a pinp document
