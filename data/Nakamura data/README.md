README
This file was generated on 2023-03-22 by Shiori Nakamura

GENERAL INFORMATION
1.	Title of Dataset: Data from: Age estimation based on blood DNA methylation levels in brown bears
2.	Author Information
Corresponding Author
Name: Dr. Hideyuki Ito3, 5 and Dr. Michito Shimozuru1
Institution: 3 Wildlife Research Center, Kyoto University, Japan; 5 Kyoto City Zoo, Japan; 1 Faculty of Veterinary Medicine, Hokkaido University, Japan
Email: 
Hideyuki Ito, itohide7@gmail.com
Michito Shimozuru, shimozuru@vetmed.hokudai.ac.jp

Co-investigators
Name: Shiori Nakamura1, Jumpei Yamazaki1, Naoya Matsumoto2, Miho Inoue-Murayama3, Huiyuan Qi3, Masami Yamanaka4, Masanao Nakanishi4, Yojiro Yanagawa1, Mariko Sashika1, and Toshio Tsubota1 
Institutions: 1 Faculty of Veterinary Medicine, Hokkaido University, Japan; 2 Noboribetsu Bear Park, Japan; 3 Wildlife Research Center, Kyoto University, Japan; 4 Shiretoko Nature Foundation, Japan; 5 kyoto City Zoo, Japan
3.	Date of data collection: 1987-2022
4.	Recommended citation for this data set: Nakamura S et al. (2023), Data from: Age estimation based on blood DNA methylation levels in brown bears, Dryad, Dataset

DATA & FILE OVERVIEW
1.	Description of data set
The data was generated to develop an age estimation model using blood DNA methylation in brown bears.

2.	File list
(1)	Name: brown_bear_blood.csv
Description: Values of the methylation levels of the samples used to build the age estimation models. This data did not include the samples that were not used for model establishment.
(2)	Name: brown_bear_blood_standardized.csv
Description: Standardized values of the methylation levels of the samples used to build the age estimation models. We used these standardized values when building the model. You need to download this file when you apply our age estimation model to a new sample.
(3)	Name: within-individual_change.csv
Description: Values of the methylation levels of the samples of wild female bears that were sampled multiple times.
(4)	Name: brown_bear_blood_test.csv
Description: Values of the methylation levels of the samples of wild female bears sampled multiple times. This data did not include the samples that were used for model establishment.
(5)	Name: brown_bear_blood_standardized_test.csv
Description: Standardized values of the methylation levels of the samples of wild female bears sampled multiple times. This data did not include the samples that were used for model establishment. This data was used to calculate MAE, Median AE, and RMSE in the "wild bears" column of Table 2 in Nakamura S et al. (2023).
(6)	Name: brown_bear_blood_all.csv
Description: Values of the methylation levels of all samples obtained in this study. This data was used to calculate MAE, Median AE, and RMSE in a SVR model without the pre-selection step.
(7)	Name: brown_bear_blood_all_standardized.csv
Description: Standardized values of the methylation levels of all samples obtained in this study. This data was used to build a SVR model without the pre-selection step.
(8)	Name: data_standardized.csv
Description: Template file for age estimation using our models. Values of methylation levels should be standardized before applying to the models. Detailed R scripts are provided in the Supplementary file_R script in Nakamura S et al. (2023).

3.	Methodological information
Methods for DNA extraction, bisulfite conversion, polymerase chain reaction and pyrosequencing were detailed in the current study (Nakamura et al. 2023) and the previous study.
・	Yamazaki, J., Meagawa, S., Jelinek, J., Yokoyama, S., Nagata, N., Yuki, M., & Takiguchi, M. (2021). Obese status is associated with accelerated DNA methylation change in peripheral blood of senior dogs. Research in Veterinary Science, 139, 193-199.

4.	DATA-specific Information
(1), (3), (4), and (6)
・	Column A (Sample_ID): Bear ID

・	Column B (birth; year, month, day): Birth date. It was assumed that all bears were born on February 1.

・	Column C (sampling_date; year, month, day): Date of the blood sampling.

・	Column D (age; year): Ages were determined at the time of blood sampling based on the assumption that all bears were born on February 1.
Reference: Friebe, A., Evans, A. L., Arnemo, J. M., Blanc, S., Brunberg, S., Fleissner, G., . . . Zedrosser, A. (2014). Factors affecting date of implantation, parturition, and den entry estimated from activity and body temperature in free-ranging brown bears. PLoS ONE, 9(7), e101410.

・	Column E (Sex): Sex, F: female, M: male.

・	Column F (environment): Growth environment (i.e., captive or wild).

・	Column G - S (CpG’s name_methylation_rate_ave; %): These columns show values of the methylation levels of the samples. As PCR for each sample was conducted in duplicate, the average value was taken as the methylation level for each sample. The CpG's names for Column G - S are SLC12A5_1, SLC12A5_2, SLC12A5_3, SLC12A5_4, POU4F2_1, POU4F2_2, POU4F2_3, POU4F2_4, VGF_1, VGF_2, VGF_3, SCGN_1, and SCGN_2, respectively.


(2), (5), and (7)
・	Column A - F are same as (1), (3), (4), and (6).

・	Column G - S (CpG’s name_methylation_rate_ave; %): These columns show standardized values of the methylation levels of the samples. The CpG's names for Column G - S are SLC12A5_1, SLC12A5_2, SLC12A5_3, SLC12A5_4, POU4F2_1, POU4F2_2, POU4F2_3, POU4F2_4, VGF_1, VGF_2, VGF_3, SCGN_1, and SCGN_2, respectively.
File (2): The values were standardized based on the data shown in the file (1).
File (5): The values were standardized based on the data shown in the file (4).
File (7): The values were standardized based on the data shown in the file (6).


(8)
・	Column A (Sample_ID): Examples of bear ID

・	Column B - N (CpG’s name_methylation_rate_ave; %): These columns show examples of standardized values of the methylation levels. The CpG's names for Column B - N are SLC12A5_1, SLC12A5_2, SLC12A5_3, SLC12A5_4, POU4F2_1, POU4F2_2, POU4F2_3, POU4F2_4, VGF_1, VGF_2, VGF_3, SCGN_1, and SCGN_2, respectively.

