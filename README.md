# MEASURING INTERNET - Guides on how to measure and analyze latency information on Internet

This repository contains two courses for students and for any person interested on learning a bit more about how to measure latency on Internet and how to analyze the data. 

The repository contains two  courses:

1. ***RIPE Atlas: Can we find where the Cloud is?***

    This course focuses on explaining what RIPE Atlas is and how we can use it to create measurements of latency on Internet.

2. ***Understanding the Cloud through data analysis***

    This course focuses on explaining how we can process and analyze the data obtained from measurements. The measurements that we use are obtained from RIPE Atlas following the first course, but the content of this course is valid for any other measurement data. 


The courses consist of three types of files: 
1. **PDF guides** that explain the details and the how-to about the Jupyter notebooks. There 
2. **Jupyter Notebooks** that contain the code that run the code that we want to generate, either to create measurements or to analyze the data. 
3. ***Data files***, which are JSON and CSV files containing the required data to run the code in the notebooks.

Note that there are two PDF guides per course: 
- One, for the student following this course, which is enough to complete the course.
- The second is aimed at the instructor of the course (if there exists an instructor guiding a group of students). This guide adds some comments on how to teach the content and few hints to do before the course happens. 

The PDF files have the title in Spanish, but the document is only in English.

## 1st course: ***RIPE Atlas: Can we find where the Cloud is?***

Please, start by reading the PDF guide file that corresponds to your role:

+ `latency1_whereIsTheCloud_GUIDE_instructor.pdf`
+ `latency1_whereIsTheCloud_GUIDE_student.pdf`

The guide will explain you the context of the courses, the motivation, and how to navigate the code. 
While reading the guide, you will be invited to open the Jupyter Notebook to follow the instructions that they contain to proceed with the course. 

There are two Jupyter Notebooks:

+ `latency1_whereIsTheCloud_A_createMeasurement.ipynb`
+ `latency1_whereIsTheCloud_A_getMeasurement.ipynb`

In order to create the measurements, you can define your own source and destination IP adresses, but we have already provided  example lists that are available at:

- `list_id_probes.csv`
- `list_url_ip_destination.csv`

## 2nd course: ***Understanding the Cloud through data analysis***

Please, start by reading the PDF guide file that corresponds to your role:

+ `latency2_cloudDataAnalysis_GUIDE_instructor.pdf`
+ `latency2_cloudDataAnalysis_GUIDE_student.pdf`

The guide will explain you the context of the courses, the motivation, and how to navigate the code. 
While reading the guide, you will be invited to open the Jupyter Notebook to follow the instructions that they contain to proceed with the course. 

There is a single Jupyter Notebook:

+ `latency2_cloudDataAnalysis.ipynb`

This notebook is self-explained and contains details that describe and comment all the steps that are applied. 

We provide the data of some RIPE Atlas measurements as an example for the course. Of course, the input data can be substituted by other measurements. The examples are provided in the JSON files: 

- `RIPE-ATLAS-measurement-nnnnnnnn.json`

where nnnnnnnn = "48819905", "48819907", "48819909" are measurements of PING packets 
and   nnnnnnnn = "48819906", "48819908", "48819910" are measurements of TRACEROUTE packets

## About

This guide has been created by the Opportunistic Architectures Lab of [IMDEA Networks Institute](https://networks.imdea.org "Developing the Science of Networks"), in the frame of the TelecoRenta Proyect from the UNICO 5G Program.        



## Contact

If you have any questions, feedback, or inquiries about the content on this repository, please feel free to reach out to our team using the following email addresses:

- For technical inquiries: [rita.ingabire@imdea.org](mailto:rita.ingabire@imdea.org)
- For technical inquiries: [antonio.bazco@imdea.org](mailto:antonio.bazco@imdea.org)


---
