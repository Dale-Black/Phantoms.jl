"""
    findCircle(point_1, point_2, point_3)

Given three points, centerpoints of the circle with [y, x] order instead of 
the more common [x, y] ordering.
"""
function find_circle(point_1, point_2, point_3)
    x1, y1 = point_1
    x2, y2 = point_2
    x3, y3 = point_3
    
    x12 = x1 - x2 
    x13 = x1 - x3  
    y12 = y1 - y2  
    y13 = y1 - y3 
    y31 = y3 - y1  
    y21 = y2 - y1
    x31 = x3 - x1  
    x21 = x2 - x1 
 
    sx13 = x1^2 - x3^2  
    sy13 = y1^2 - y3^2
    sx21 = x2^2 - x1^2  
    sy21 = y2^2 - y1^2  
  
    f = (((sx13) * (x12) + (sy13) * (x12) + (sx21) * (x13) + (sy21) * (x13)) ÷ (2 * ((y31) * (x12) - (y21) * (x13)))) 
              
    g = (((sx13) * (y12) + (sy13) * (y12) + (sx21) * (y13) + (sy21) * (y13)) ÷ (2 * ((x31) * (y12) - (x21) * (y13))))  
  
    # eqn of circle be x^2 + y^2 + 2*g*x + 2*f*y + c = 0 where center is (h = -g, k = -f)  
    center_insert = [-g, -f]

    return center_insert
end