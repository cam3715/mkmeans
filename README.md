# mkmeans

Copyright (C) 2016 Camden Bock, Bates College.

This program is modified from Ahmad Alsahaf`s package: amjams/mixedkmeans (Licensed under GNU GPL).
This program is licensed under GNU GPL v3.0.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.\\

The source code is freely available at:\\ \url{https://github.com/cam3715/mkmeans}\\


Cite this source code with the following paper:\\

\textbf{MLA:}\\
Bock, Camden Glenn, and Bonnie J. Shulman. "Mixed K-means Clustering in Computer Adaptive Learning." \emph{SCARAB} (2016) \emph{Department of Mathematics,} Bates College, May 2016. Web. \href{url}{http://scarab.bates.edu/}.\\

\textbf{APA:}\\
Bock, C. G., \& Shulman, B. J. (2016). Mixed K-means clustering in computer adaptive learning. \emph{SCARAB}. Retrieved from \href{url}{http://scarab.bates.edu/}.\\

You should have received a copy of the GNU General Public License
along with this program.  If not, see \href{url}{http://www.gnu.org/licenses/}.

##Abstract
The ASSISTments project from Worcester Polytechnic Institute provides a free web-based intelligent tutoring system including two levels of differentiation, that are manually programmed by teachers and 
researchers.
Problems assigned through ASSISTments can be programmed in trees, where the sequence of problems adapts to the student's performance on each question.
Within each problem, if a student enters an incorrect response the ASSISTments system provides scaffolded feedback to target the student's misconception.
This thesis begins to develop an educational data mining algorithm to automate this differentiation.
First, an adaption of Alsahaf's mixed k-means clustering algorithm is proposed to handle a mix of categorical and numeric data.
Second, the algorithm is implemented in MATLAB and its performance is compared to Alsahaf's results on benchmark data sets.
Finally, the MATLAB implementation is applied to ASSISTments data sets from 2009 and 2012 to develop a predictive model.


##Introduction
Intelligent tutoring systems provide instruction that is differentiated to individual student needs.  
Computer adaptive learning systems, a subset of intelligent tutoring systems, differentiate for individual students by frequently adapting to a student's behavior.
While other intelligent tutoring systems may rely on tracking, defining students by demographic labels (e.g. socioeconomic status, ethnicity, IEP or special education designation), computer adaptive learning has the potential for more equitable instruction.
One of the motivations for the work in this thesis is to help realize that potential.

Intelligent tutoring systems also have the potential to provide a large number of students with additional support structures and more equitable access to content.
Standards based education (SBE) requires frequent and specific intervention and allows for multiple methods for demonstrating proficiency.
As a support to the teacher and student, intelligent learning systems can provide struggling students with continuously differentiated skills practice.  
Additionally, proficiency on standards is assessed while the student is simultaneously learning.
This eases the burden on teachers and schools of reassessing failing students both by minimizing time spent grading and generating alternative assessments.

Finally, many intelligent learning systems are web-based, and can be accessed by netbooks, smartphones, tablets and chromebooks in addition to traditional desktops and laptops.
With 1-1 technology initiatives, these intelligent learning systems can be accessed by students without additional hardware or software costs to the school districts.

ASSISTments, a web-based intelligent tutoring system from Worcester Polytechnic Institute (WPI), funded by the National Science Foundation, is a free resource for districts, educators and researchers.  
Although ASSISTments has adaptive differentiation, it requires manual entry by the user.
Commercial systems have automated this differentiation of problem sets.
The addition of computer adaptive features to the ASSISTments platform would increase its value to educators and researchers and decrease the labor cost of implementing ASSISTments in local districts.
Computer adaptive features would be most useful in either the design of differentiated and adaptive sequences of problems towards proficiency on a skill or in diagnosing misconceptions and appropriate interventions in students' responses.

ASSISTments has collected data from thousands of students over the last decade, focused on secondary mathematics in Maine and Massachusetts.
\footnote{A recent study by WPI, SRI and the University of Maine note both significant gains in the performance of grade seven students in a randomized control trial including 44 schools in Maine.}
ASSISTments currently allows teachers to assign linear sequences of problems, or a sequence of problems determined by the student's performance (correct/incorrect) on the previous problem.
Creating these problem sets is a time-intensive task for the teacher, and does not take advantage of years of student data.
ASSISTments provides scaffolding, an additional layer of differentiation  within each problem.
When a student enters an incorrect response, the system diagnoses the student's error, and asks additional questions that break the problem into smaller pieces around that student's error.
If the student has sufficient background knowledge/instruction, this approach should identify the student's misconceptions and correct them with additional questions within the students zone of proximal development.

In order to automate either the problem sequence or problem scaffolding, students or their responses need to be clustered, so that regression models can be made on the smaller partitions, for more precise predictions.
Clustering uses unsupervised machine learning, where the k-means algorithm efficiently partitions data that is structured in $m-$dimensional spheres.
Points near a cluster center (mean) can be predicted to behave as the cluster center or based on the deviation from the cluster center.

To find spherical groupings, k-means clustering needs a well-defined distance between data points.
The data from ASSISTments, as well as many other systems, has a mix of categorical and numeric data.
The squared Euclidean distance cannot be directly applied to categorical data.
Ahmad et al. propose a measure of distance between categorical variables, and a weighted measure of distance between data points with numeric and categorical attributes \cite{ahmad_k-mean_2007}.
Alsahaf implemented Ahmad's algorithm into MATLAB, tailored to a specific set of data \cite{alsahaf_amjams_2015}.

In this thesis, we generalize Alsahaf's mixed k-means algorithm.  
Our algorithm accepts dense matrices as an input.
The results of our algorithm agree with those reached by Ahmad on benchmark data sets.
