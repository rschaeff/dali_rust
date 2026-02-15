"""
WOLF module: SSE-based protein structure comparison.

Reimplementation of the WOLF stage from DaliLite.v5.
"""

from .wolf import wolf_compare, WolfResult

__all__ = ['wolf_compare', 'WolfResult']
