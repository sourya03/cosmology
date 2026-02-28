from flask import Flask, request, jsonify, render_template
from flasklogic import calculate_cosmology   # your earlier Python logic

app = Flask(__name__)

@app.route('/')
def index():
    return render_template('index.html')  # the pretty frontend HTML

@app.route('/cosmology', methods=['POST'])
def cosmology_api():
    data = request.json
    H0 = data['H0']
    Omega_m = data['Omega_m']
    Omega_lambda = data['Omega_lambda']
    z = data['z']
    
    result = calculate_cosmology(H0, Omega_m, Omega_lambda, z)
    return jsonify(result)

if __name__ == '__main__':
    app.run(debug=True)
