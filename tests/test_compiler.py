import unittest
from subprocess import check_output


class TestCompiler(unittest.TestCase):
    def test_cmake(self):
        if "no g++" in check_output(["which", "g++"]).decode():
            raise EnvironmentError("cmake not found")

if __name__ == '__main__':
    unittest.main()
