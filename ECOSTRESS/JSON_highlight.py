"""
This module handles syntax highlighting for console output. This makes console notifications much more readable.

Developed by Gregory Halverson in the Jet Propulsion Laboratory Year-Round Internship Program
(Columbus Technologies and Services), in coordination with the ECOSTRESS mission
and master's thesis studies at California State University, Northridge.
"""

import json
from builtins import str

from pygments import highlight
from pygments.formatters.terminal import TerminalFormatter
from pygments.lexers.data import JsonLexer

__author__ = 'Gregory Halverson'


def JSON_highlight(tree: dict, indent: int = 2) -> str:
    """
    Converts python dictionary tree mapping to JSON and highlights JSON syntax for console output.
    :param tree: dictionary tree/mapping
    :param indent: number of spaces to use for indents
    :return: JSON-highlighted string
    """
    if type(tree) is dict:
        tree = json.dumps(tree, indent=indent)

    tree = highlight(str(tree), JsonLexer(), TerminalFormatter())

    return tree
