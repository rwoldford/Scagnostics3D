Instructions for testing scagnostics3D
Aug 13, 2009

NOTE: 
=====
This package includes Qhull files, version 2003.1 for Windows, 2003.1 for Unix (see http://www.qhull.org/) for the purpose of test. Please refer to Qhull official website for any documentations.

The files in folder \Mac are for Macintosh, \Windows are for Windows XP.  


Prerequisites:
=============
R2.7.1 or R2.8.0 is recommended. Please make sure you have rJava package installed. Other useful packages might be rggobi and scatterplot3d. 


List of files:
=============
README

scag.R
testdata.R

\Mac\scagnostics3D.jar
Mac\qhull\cmdfile
Mac\qhull\qdelaunay
Mac\qhull\qconvex
Mac\qhull\.lib\*.*  (hidden files, don NOT change them)
Mac\qhull\.deps\*.*  (hidden files, don NOT change them)

\Windows\qdelaunay.exe
\Windows\qconvex.exe
\Windows\scagnostics3D.jar
==============


Steps to test scagnostics3D on Mac:
===================================
1. Create folder /tmp/scagnostics3D/ (NOTE: must be this folder, cannot be changed)

2. Copy /qhull(directory) and scagnostics3D.jar to the folder created in step 1

3. You can put scag.R and testdata.R anywhere you are used to. 

3. Launch R and open scag.R, run the entire script.

4. At the end of running, if you read some numbers showing the scagnostics info(from R console), congratulations, you are done with the setting. 

5. Now you can try testdata.R, or other data set as long as it meets the scagnostcis3D required input format.


Steps to test scagnostics3D on Windows XP:
=========================================

1. Create folder C:\scagnostics3D  (or you can create any folder you like)

2. Copy all the files in \Windows to the folder created in step 1.

3. scag.R and testdata.R can be put any directory.

4. Launch R and open scag.R, modify the script to uncomment the line #94:
 .jaddClassPath("C:/scagnostics3D/scagnostics3D.jar")
Change the directory as necessary, if you create one other than recommended in step 1.
Of course, you need to comment the line for Mac (#91). 

5. Run the entire script until a bunch of scagnostics numbers appear in the R console. 

