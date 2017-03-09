from myreport import html_report
import random

# PUT ALL FUNCTIONS HERE
#def draw_ball():

#def draw_ball2():


# MAIN

myreport=html_report("s6_task2.html")
myreport.add_heading("2 Generating a random choice")

# A
myreport.add_subheading ("a.")
myreport.add_text("")

# B
myreport.add_subheading ("b.")
myreport.add_text("green was  %")
myreport.add_text("red was  %")

# C
myreport.add_subheading ("c.")
myreport.add_text("Result of function draw_ball")
#result=draw_ball()
#myreport.add_text(result)

# D
myreport.add_subheading ("d.")
myreport.add_text("Result of 10,000 draws")

# E
myreport.add_subheading ("e.")
myreport.add_text("Result of function draw_ball2(10,3)")
#result=draw_ball2(10,3)
#myreport.add_text(result)




# ADD SOURCE CODE AND WRITE REPORT
# THEN OPEN IT IN A BROWSER

# INSERT SOURCE CODE INTO REPORT
myreport.add_subheading('Python Code')
myreport.add_source(__file__)

# WRITE REPORT TO FILE

myreport.write()


# OPEN FILE IN BROWSER

myreport.view()
