# Lam_Morsut_GJSM

These are the original run codes for the simulations in: https://doi.org/10.1021/acssynbio.0c00369. 
They are numbered according to the figure label. (i.e. F7B=Figure 7B).

Getting Started
1) First get Compucell3D v.3.7.8 (CC3D) from the CompuCell3D website. 
Here is the link for the installation files, but this may change over time: https://sourceforge.net/projects/cc3d/files/3.7.8/
2) Install CompuCell3D as per the installation file.

Running the code
1) Download the codes in this repository. You can download them individually or use the code tab.
![image](https://user-images.githubusercontent.com/68087210/165001891-0da7a892-3d12-49b9-ae3c-3d435f489eeb.png)
2) Once you have the codes downloaded and CC3D installed, you can run either Twedit++ or CompuCell3D player directly.
My installation is in windows. I've included my directory in the image but your may vary.
![image](https://user-images.githubusercontent.com/68087210/165002024-67410987-7c05-4f71-9cd2-c9d52e110a0a.png)

Running from CC3D player. 
1) Click to run the CompuCell3D file, then click file on the top left, and open simulation file.

![image](https://user-images.githubusercontent.com/68087210/165002116-d0a2b659-af77-4760-9977-615f0aabdd57.png)

2) Navigate to where your simulation file is. Then click open.

![image](https://user-images.githubusercontent.com/68087210/165002199-dde1354d-6515-410f-9651-950e290a7c09.png)

3) You'll enter the next level of the folder. Select the CC3D file and open it. 

 ![image](https://user-images.githubusercontent.com/68087210/165002253-707d6191-d399-42c1-840b-31a4d2bff0a3.png)
 
4) The simulation file is loaded. Next you can run the file by clicking the play button.

![image](https://user-images.githubusercontent.com/68087210/165002335-caeac435-feb3-4468-89dc-e61eecb244d7.png)

5) The code will output files. To find where this output is: go to Tools in the player, then configuration.

![image](https://user-images.githubusercontent.com/68087210/165002471-2a92e6a0-1fd1-4d18-8117-96feb111e981.png)

6) Then go to the output tab, and look at the directory for output. You can also change to where you want the output. Mine is set to my desktop. 
Don't forget to hit Apply so it sticks.

![image](https://user-images.githubusercontent.com/68087210/165002564-eddf3425-487a-46d7-880f-0b895f417a75.png)

The player should run and generate your output files in your wanted directory :)

Running from Twedit++. 
1) In step two of "Running from CC3D player" open twedit++ instead.
2) Open your file by going to CC3D project on the top left, then select Open CC3D Project

![image](https://user-images.githubusercontent.com/68087210/165002718-4016b1d5-a8d0-47a8-9ace-326e55bb7c81.png)

3) Open your project.

![image](https://user-images.githubusercontent.com/68087210/165002777-60aa6a1c-13a3-4f61-b9f2-a1336a100862.png)
![image](https://user-images.githubusercontent.com/68087210/165002831-12bb31eb-9f33-4636-8f5e-27eb702c9fc5.png)

4)Double click the simulation name. This will expand the code to its separate file parts.

![image](https://user-images.githubusercontent.com/68087210/165002919-f18a4a31-a0a1-4c94-becb-e63dc85ee40c.png)

5) You can run the code by right clicking the simulation file, then "Open in Player"

![image](https://user-images.githubusercontent.com/68087210/165003001-20d54de9-2532-41e2-b8e8-906fd3c4f412.png)

Troubleshooting.
If the simulation crashes shortly after running, it may be becasue the simulation is set to 8 cores for parallelization, but your system cannot support that.
1) In this case, from the twedit instruction set (From Step 4), navigate to the tab 'ELUGMSteppables'.
![image](https://user-images.githubusercontent.com/68087210/165003098-c78b42ed-4874-401e-a68d-989ccc51c281.png)
2) Go to line 45, set USEDNODES=1 (Or whatever is appropriate for your sysem)
) ![image](https://user-images.githubusercontent.com/68087210/165003135-19743244-4bc2-48bb-803a-5c7a02446d25.png)
I find 8 offers a good balance of performance and speed, without splitting the lattice into too many pieces.

Contact me anytime for help or with questions calvin.lam.k@gmail.com

