from flask import Flask

app = Flask(__name__) # create an instance of the Flask class

def index():
    return "中国联通"