from flask import Flask, render_template, request
import subprocess
import sys
# import main

app = Flask(__name__)

print("For best results, please use Google Chrome for browser-based GUI")


@app.route('/')
def index():
    return render_template('index.html')

@app.route('/result', methods=['POST', 'GET'])
def get_result():
    if request.method == 'POST':
        result = request.form
        data = result.to_dict()
        # print(data)
        f = open('./userdata.txt', 'w')
        for p in data:
            f.write(f'{data[p]}\n')
        f.close()

        subprocess.run([sys.executable, "main.py"])
        
        # change the second parameter below to the data we analyse in main.py.
        # after the main computation we will display the useful output data on result.html using templating
        return render_template("result.html", result=result)



if __name__ == '__main__':
    app.run(port=3000, debug=True)
