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
