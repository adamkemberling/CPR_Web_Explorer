Continuous Plankton Recorder Web Explorer App
================
Adam A. Kemberling
9/27/2020

# About the App

This is repository is for the development and deployment of an R Shiny
web application for the purpose of visualizing and exploring continuous
plankton recorder (CPR) survey data.

To build and run the docker image:

> docker build –tag gmri/bashful-badger:1.0 .  
> docker run -p 3838:3838 –name cpr gmri/bashful-badger:1.0
