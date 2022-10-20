from flask import Flask, render_template, request
app = Flask(__name__)

posts = [
    {
        'author': 'Corey Shafer',
        'title': 'Blog Post 1',
        'content': 'First post!',
        'date_posted': 'April 20, 2018',
    },
    {
        'author': 'Jane Doe',
        'title': 'Blog Post 2',
        'content': 'Dogs r Cute',
        'date_posted': 'April 21, 2018',
    }
]


@app.route("/")
@app.route("/home")
def home():
    return render_template("index.html", posts=posts)


@app.route("/about")
def about():
    return "<h1>About Page</h1>"


if __name__ == '__main__':
    app.run(debug=True)
