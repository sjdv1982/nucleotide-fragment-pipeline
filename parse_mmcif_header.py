import itertools
import re

parenthesis_match = re.compile(r"\([^\(\)]*?\)")


def _parseOperationSubExpression(expression) -> list[int]:
    expression = expression.strip()
    operations = []
    if expression.find(",") > -1:
        # comma-separated list of expressions
        terms = expression.split(",")
        for term in terms:
            operations += _parseOperationSubExpression(term)
    elif expression.find("-") > -1:
        # range
        start, end = expression.split("-")
        start = int(start)
        end = int(end)
        for n in range(start, end + 1):
            operations.append(str(n))
    else:
        # single operation
        operations = [expression]
    return operations


def _parseOperationExpression(expression):
    subexpressions = [m[1:-1] for m in re.findall(parenthesis_match, expression)]
    if not subexpressions:
        subexpressions = [expression]
    parsed_subexpressions = [
        _parseOperationSubExpression(subexpression) for subexpression in subexpressions
    ]
    unrolled_subexpressions = list(itertools.product(*reversed(parsed_subexpressions)))
    return unrolled_subexpressions


def _nest_dict(cif):
    result = {}
    for k, v in cif.items():
        parent = result
        fields = k.split(".")
        for field in fields[:-1]:
            if field not in result:
                result[field] = {}
            parent = result[field]
        parent[fields[-1]] = v
    return result


def _replace_operation_expressions(d):
    for k, v in d.items():
        if isinstance(v, dict):
            _replace_operation_expressions(v)
            continue
        if k == "oper_expression":
            vlist = v
            if isinstance(v, str):
                vlist = [v]
            vlist_parsed = []
            for vv in vlist:
                vv_parsed = _parseOperationExpression(vv)
                vlist_parsed.append(vv_parsed)
            d[k] = vlist_parsed


def parse_mmcif_header(mmcif_data):
    from io import StringIO
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict

    mmcif_obj = StringIO(mmcif_data)
    cif = MMCIF2Dict(mmcif_obj)
    result = {}
    if "data_" in cif:
        result["data"] = cif["data_"]
    for k in cif:
        if k == "data_":
            continue
        kk = k[1:]
        if k[0] != "_":  # should not happen
            kk = k
        if kk.startswith("atom_site"):
            continue
        result[kk] = cif[k]
    result_nested = _nest_dict(result)
    _replace_operation_expressions(result_nested)
    return result_nested


if __name__ == "__main__":
    import sys
    import json

    mmcif_file = sys.argv[1]
    outfile = sys.argv[2]
    data = parse_mmcif_header(open(mmcif_file).read())
    with open(outfile, "w") as f:
        json.dump(data, f)
