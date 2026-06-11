import json

with open("madspace.json") as f:
    data = json.load(f)

color_flows = {}
color_orders = {}
helicities = {}
pdg_ids = {}

def dict_index(d, k):
    if k in d:
        return d[k]
    else:
        index = len(d)
        d[k] = index
        return index

for chan in data["channels"]:
    del chan["channel"]
    chan["phasespace_order"] = [c - 1 for c in chan["phasespace_order"]]
    for proc in chan["processes"]:
        proc["color_flows"] = dict_index(
            color_flows, tuple(zip(proc["color_flows1"], proc["color_flows2"]))
        )
        del proc["color_flows1"]
        del proc["color_flows2"]
        del proc["process"]
        proc["color_order"] = dict_index(color_orders, tuple(c - 1 for c in proc["color_order"]))
        proc["multichannels"] = [mc - 1 for mc in proc["multichannels"]]
        for mat in proc["matrix_elements"]:
            mat["pdg_ids"] = dict_index(pdg_ids, tuple(mat["pdg_ids"]))
            del mat["label"]
        for hel_group in proc["helicities"]:
            hel_group[:] = [
                dict_index(helicities, tuple(hel)) for hel in hel_group
            ]

out_dict = {
    "channels": data["channels"],
    "color_flows": list(color_flows.keys()),
    "color_orders": list(color_orders.keys()),
    "pdg_ids": list(pdg_ids.keys()),
    "helicities": list(helicities.keys()),
    "xml_header": data["xml_header"], 
}
with open("madspace_converted.json", "w") as f:
    json.dump(out_dict, f)
