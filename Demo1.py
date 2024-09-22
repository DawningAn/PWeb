from flask import Flask, render_template

app = Flask(__name__)  # create an instance of Flask class


# 创建了一个网址 /show/info 和一个函数名 index() 的对应关系
# 以后用户在浏览器中访问 /show/info 这个网址时，会自动调用 index() 函数
@app.route("/show/info")
def index():
    # return "忘川如斯"
    # Flask 内置了一个模板引擎，内部会自动加载 templates 文件夹下的文件并读取内容，将内容返回给用户
    # 默认会去当前项目的 templates 文件夹下寻找模板文件
    return render_template("index.html")  # 调用模板文件


@app.route("/qs/list")
def qs_list():
    return render_template("qs_list.html")

@app.route("/zzbj/list")
def zzbj_list():
    return render_template("zzbj_list.html")

@app.route("/register")
def register():
    return render_template("register.html")

if __name__ == '__main__':
    app.run()  # 启动 Flask 服务器
