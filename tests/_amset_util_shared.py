"""测试用：amset.util 的共享数组三个函数（与 AMSET 0.4.19 / 0.5.1 逐字相同，见 amset/util.py）。"""
from multiprocessing.sharedctypes import RawArray

import numpy as np


def create_shared_array(data: np.ndarray, return_shared_data=False):
    data = np.asarray(data)
    if data.dtype == np.complex128:
        data_type = "complex"
        data_buffer = RawArray("d", int(np.prod(data.shape)) * 2)
    else:
        data_type = np.ctypeslib.as_ctypes_type(data.dtype)
        data_buffer = RawArray(data_type, int(np.prod(data.shape)))

    buffer = (data_buffer, data.shape, data_type)
    data_shared = array_from_buffer(buffer)
    data_shared[:] = data[:]

    if return_shared_data:
        return buffer, data_shared
    else:
        return buffer


def array_from_buffer(buffer):
    data_buffer, data_shape, data_type = buffer
    if data_type == "complex":
        return np.frombuffer(data_buffer).view(np.complex128).reshape(data_shape)
    else:
        return np.frombuffer(data_buffer, dtype=data_type).reshape(data_shape)


