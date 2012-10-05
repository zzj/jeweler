import os

os.system("echo config week 8 4 > command_list")
os.system("bash test.sh >> command_list")
os.system("echo pause >> command_list")
os.system("echo config week 8 1 >> command_list")
os.system("python3 combine.py >> command_list")
os.system("echo pause >> command_list")
os.system("echo config week 8 1 >> command_list")
os.system("python3 merge.py >> command_list")




