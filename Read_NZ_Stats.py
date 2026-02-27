import struct

def Read_NZ_Stats (filename):
    with open(filename, "rb") as f:
        # ---------- Read header ----------
        header_fmt = "<iiIIfffdQ"
        header_size = struct.calcsize(header_fmt)
    
        header_data = f.read(header_size)
        L, NQ, state_i, state_f, aR, aI, X_magic, total_paths, total_NZ_paths = struct.unpack(header_fmt, header_data)

        header = {}
        header["L"] = L
        header["NQ"] = NQ
        header["state_i"] = state_i
        header["state_f"] = state_f
        header["aR"] = aR
        header["aI"] = aI
        header["X_magic"] = X_magic
        header["total_paths"] = total_paths
        header["total_NZ_paths"] = total_NZ_paths

        # ---------- Prepare record format ----------
        record_fmt = "<" + ("I" * L) + "ff"
        record_size = struct.calcsize(record_fmt)

        #print("Each record size:", record_size, "bytes")

        # ---------- Read records ----------
        paths = []
        pathsW = []

        for _ in range(total_NZ_paths):
            data = f.read(record_size)
            if len(data) != record_size:
                raise ValueError("Unexpected end of file")

            unpacked = struct.unpack(record_fmt, data)

            uints = unpacked[:L]
            floats = unpacked[L:]

            paths.append(uints)
            pathsW.append(floats)
    return header, paths, pathsW