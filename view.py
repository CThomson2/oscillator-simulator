from flask import Flask, render_template, request
# import main

app = Flask(__name__)

print("For best results, please use Google Chrome for browser-based GUI")


@app.route('/')
def student():
    return render_template('index.html')


@app.route('/result', methods=['POST', 'GET'])
def result():
    if request.method == 'POST':
        result = request.form
        data = result.to_dict()
        for param in data:
            print(data[param])
        # main.main_function(data)
        return render_template("result.html", result=result)


if __name__ == '__main__':
    app.run(debug=True)
