import os
from pathlib import Path

def check_path_like(obj):
    try:
        # This will raise a TypeError if obj is not path-like
        path_str = os.fspath(obj)
        return f"This is a path-like object: {path_str}"
    except TypeError:
        return "This is not a path-like object."

# Example usages
path_obj = Path('path/to/file.txt')
string_obj = 'Just a normal string'
bytes_buffer = b'Some binary data'

print(check_path_like(path_obj))  # Output: This is a path-like object.
print(check_path_like(string_obj)) # Output: This is a path-like object.
print(check_path_like(bytes_buffer)) # Output: This is not a path-like object.
