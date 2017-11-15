import unittest
import math
import io 
import os

import precision_recall

class TestComputePrecisionRecall(unittest.TestCase):

    
    def test_compute_precision_recall_point(self):
        retrieved = set(["a", "b", "c"])
        relevant = set(["a", "b", "d", "f"])

        point = precision_recall.compute_precision_recall_point(
            retrieved, relevant)

        expected = [.6666, .5]

        self.assertTrue(
            math.isclose(point[0], expected[0], abs_tol=0.0001))
        self.assertTrue(
            math.isclose(point[1], expected[1], abs_tol=0.0001))
    

    def test_compute_precision_recall_curve(self):  
        ranked_items = [set("a"),
                        set(["b", "c"]), 
                        set(["d", "e"])]

        relevant_items = set(["a","b","e"])

        curve_points = precision_recall.compute_precision_recall_curve(
            ranked_items, relevant_items)

        expected = [(1, 0),
                    (1, .3333),
                    (.6666, .6666),
                    (.6, 1)]
        
        for pair in zip(curve_points, expected):
            self.assertTrue(
                math.isclose(pair[0][0], pair[1][0], abs_tol=0.0001))
            self.assertTrue(
                math.isclose(pair[0][1], pair[1][1], abs_tol=0.0001))


class TestPrecisionRecallIO(unittest.TestCase):
    
    def test_write_precision_recall_values(self):
        points = [(1,0), (1, .3333), (.6666, .6666), (.6, 1)]

        output = io.StringIO()
        precision_recall.write_precision_recall_values(output, points)

        expected = "1,0%s1,0.3333%s0.6666,0.6666%s0.6,1%s" % \
            (os.linesep, os.linesep, os.linesep, os.linesep)

        self.assertEqual(expected, output.getvalue())


    def test_read_precision_recall_values(self):
        infile = io.StringIO("1,0%s1,0.3333%s0.6666,0.6666%s0.6,1%s" % 
            (os.linesep, os.linesep, os.linesep, os.linesep))

        expected = [(1,0), (1, .3333), (.6666, .6666), (.6, 1)]

        points = precision_recall.read_precision_recall_values(infile)

        self.assertEqual(points, expected)


if __name__ == '__main__':
    unittest.main()
