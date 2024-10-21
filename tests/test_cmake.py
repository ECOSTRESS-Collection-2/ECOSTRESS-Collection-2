import unittest
from subprocess import check_output


class Testcmake(unittest.TestCase):
    def test_cmake(self):
        if "no cmake" in check_output(["which", "cmake"]).decode():
            raise EnvironmentError("cmake not found")

if __name__ == '__main__':
    unittest.main()
