| Date     | Description |
| :---        | :---   |
| 2021-04-16  | created and set up the repo of the project on Kyiv urban heat island (UHI) long-term spatio-temporal transformations |
| 2021-04-19  | downloaded some raw data for testing |
| 2021-04-21  | learned to unpack binaries |
| 2021-05-03  | decided on L8 binaries interpretation  to include only single bit values decoded by bits 1:6, excluding clear and water pixels (bits 7:8), and excluding confidence level values decoded by bits 9:16 |
| 2021-05-04  | learned to create binary masks based on decoded L8 QA values to exclude invalid pixels from further analysis |
| 2021-05-05  | some code clean-up  |
| 2021-05-06  | testing and prepping the code for automated processing |
| 2021-05-07  | organized code into a function: some files are not properly processed - needs investigation |
| 2021-05-12  | since stuck with cropping, proceeded so far to the climate data records analysis |
| 2021-05-13 | analysis of climate data records to find consecutive NAs that may affect the results |
| 2021-05-14 | developed a workflow to analyze growing season parameters excluding years with March-October having 5 and more consecutive days of missing observations  |
| 2021-05-23 | downloaded 1881-2020 data, introduced linear interpolation via zoo library to fill in data gaps less than 5 days |