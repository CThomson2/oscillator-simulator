import math
import numpy as np
import matplotlib.pyplot as plt

class Methods():
    
    # Root finding methods
    def newton():
        pass

    def muller(f, a, b, c, k, max_iter): 
  
        res = 0; 
        i = 0; 
    
        while (True): 
        
            # Calculating various constants required to calculate x3 in
            # Muller's Formula
            f1 = f(a, k); f2 = f(b, k); f3 = f(c, k); 
            d1 = f1 - f3;  
            d2 = f2 - f3; 
            h1 = a - c;  
            h2 = b - c; 
            a0 = f3; 
            a1 = (((d2 * pow(h1, 2)) - 
                (d1 * pow(h2, 2))) / 
                ((h1 * h2) * (h1 - h2))); 
            a2 = (((d1 * h2) - (d2 * h1)) / 
                ((h1 * h2) * (h1 - h2))); 
            x = ((-2 * a0) / (a1 + 
                abs(math.sqrt(a1 * a1 - 4 * a0 * a2)))); 
            y = ((-2 * a0) / (a1 - 
                abs(math.sqrt(a1 * a1 - 4 * a0 * a2)))); 
    
            # Taking the root which is closer to x2 
            if (x >= y): 
                res = x + c; 
            else: 
                res = y + c; 
            # checking for resemblance of x3 with x2 till two decimal places 
            m = res * 100; 
            n = c * 100; 
            m = math.floor(m); 
            n = math.floor(n); 
            if (m == n): 
                break; 
            a = b; 
            b = c; 
            c = res; 
            if (i > max_iter): 
                print("Root cannot be found using",  
                                "Muller's method"); 
                break; 
            i += 1; 
        if (i <= max_iter): 
            # print("The value of the root is", 
            #                 round(res, 4)); 
            # print(i)
            return res

    # ODE solvers
    def euler_simple():
        return 2 + 2

    def runge_kutta():
        pass

    def implicit_euler():
        pass

# print(Methods.euler_simple())