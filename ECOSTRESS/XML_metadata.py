def write_XML_metadata(metadata_dict: dict, filename: str):
    XML_string = ""
    XML_string += '<?xml version="1.0" encoding="UTF-8"?>\n'
    XML_string += '<cas:metadata xmlns:cas="http://oodt.jpl.nasa.gov/1.0/cas">\n'

    for key, value in metadata_dict.items():
        XML_string += '   <keyval type="vector">\n'
        XML_string += f"      <key>{key}</key>\n"
        XML_string += f"      <val>{value}</val>\n"
        XML_string += "   </keyval>\n"

    XML_string += "</cas:metadata>\n"

    with open(filename, "w") as file:
        file.write(XML_string)
