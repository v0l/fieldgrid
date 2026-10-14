import numpy as np
import openseespy.opensees as ops

E, NU = 30e9, 0.2
G = E / (2 * (1 + NU))
H, BX, BY, STOREYS = 3.2, 6.0, 5.0, 3
FLOOR_MASS = 20000.0


def rect(w, d):
    a, b = max(w, d), min(w, d)
    j = a * b**3 * (1 / 3 - 0.21 * b / a * (1 - b**4 / (12 * a**4)))
    return w * d, d * w**3 / 12, w * d**3 / 12, j


COLUMN, BEAM_X, BEAM_Y = rect(0.4, 0.4), rect(0.3, 0.5), rect(0.25, 0.45)


def build():
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    tags = {}
    for k in range(STOREYS + 1):
        for j in range(3):
            for i in range(3):
                tag = len(tags) + 1
                tags[(i, j, k)] = tag
                ops.node(tag, i * BX, j * BY, k * H)
                if k == 0:
                    ops.fix(tag, 1, 1, 1, 1, 1, 1)
                else:
                    ops.mass(tag, FLOOR_MASS, FLOOR_MASS, FLOOR_MASS, 0, 0, 0)
    members = []
    for k in range(STOREYS):
        for j in range(3):
            for i in range(3):
                members.append((tags[(i, j, k)], tags[(i, j, k + 1)], COLUMN, (1, 0, 0)))
    for k in range(1, STOREYS + 1):
        for j in range(3):
            for i in range(2):
                members.append((tags[(i, j, k)], tags[(i + 1, j, k)], BEAM_X, (0, 0, 1)))
        for j in range(2):
            for i in range(3):
                members.append((tags[(i, j, k)], tags[(i, j + 1, k)], BEAM_Y, (0, 0, 1)))
    for e, (a, b, (area, iy, iz, j), up) in enumerate(members, 1):
        pa, pb = np.array(ops.nodeCoord(a)), np.array(ops.nodeCoord(b))
        x = (pb - pa) / np.linalg.norm(pb - pa)
        u = np.array(up, float)
        y = u - x * u.dot(x)
        y /= np.linalg.norm(y)
        z = np.cross(x, y)
        ops.geomTransf("Linear", e, *z)
        ops.element("elasticBeamColumn", e, a, b, area, E, G, j, iy, iz, e)
    return tags, members


tags, members = build()
ops.timeSeries("Linear", 1)
ops.pattern("Plain", 1, 1)
for k in range(1, STOREYS + 1):
    ops.load(tags[(0, 0, k)], 10000.0 * k, 0, 0, 0, 0, 0)
ops.system("BandGeneral")
ops.numberer("RCM")
ops.constraints("Plain")
ops.integrator("LoadControl", 1.0)
ops.algorithm("Linear")
ops.analysis("Static")
ops.analyze(1)
print("roof", [ops.nodeDisp(tags[(i, j, STOREYS)], 1) for (i, j) in [(0, 0), (2, 2)]])
print("roof_twist", ops.nodeDisp(tags[(1, 1, STOREYS)], 6))
ops.reactions()
print("base_moment", ops.nodeReaction(tags[(0, 0, 0)], 5))
ops.wipeAnalysis()
build()
values = ops.eigen("-fullGenLapack", 6)
print("periods", [2 * np.pi / np.sqrt(v) for v in values])
