import math
import os
import subprocess

if __name__ == "__main__":
    
    c = 500
    r = 50.0
    x = 0.0
    y = 0.0
    
    # create edge file
    with open("my.edg.xml", "w") as output:
 
        angle = 2 * math.pi / c
        shape1 = ["%.2f,%.2f" % (math.cos(i * angle) * r + x,
                                math.sin(i * angle) * r + y) for i in range(round(c/2)+1)]
        
        shape2 = ["%.2f,%.2f" % (math.cos(i * angle) * r + x,
                                math.sin(i * angle) * r + y) for i in range(round(c/2), c+1)]
        print('''
    <edges>
        <edge id="edge1" from="n1" to="n2" shape="%s" numLanes="2" speed="20"/>
        <edge id="edge2" from="n2" to="n1" shape="%s" numLanes="2" speed="20"/>
    </edges>''' % (" ".join(shape1), " ".join(shape2)), file=output)
    
    # create node file    
    with open("my.nod.xml", "w") as output:
 
        x1 = "%.2f" % r
        y1 = "%.2f" % 0.0
        x2 = "%.2f" % -r
        y2 = "%.2f" % 0.0
        print('''
    <edges>
        <node id="n1" x="%s" y="%s" type="priority"/>
        <node id="n2" x="%s" y="%s" type="priority"/>
    </edges>''' % (x1, y1, x2, y2), file=output)
    
    # create the net file
    os.system("netconvert --node-files=my.nod.xml --edge-files=my.edg.xml --output-file=CircularRoad.net.xml")
    
    
    