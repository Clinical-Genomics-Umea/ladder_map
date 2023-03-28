import pandas as pd
import pandera as pa
import yaml
from yaml.loader import SafeLoader

with open('fragment_analyzer/config/method_validation_criteria.yaml') as f:
    criteria = yaml.load(f, Loader=SafeLoader)

schema = pa.DataFrameSchema({
    "name": pa.Column(str,
                      checks=pa.Check.str_matches(criteria['regex_name']),
                      required=True,
                      unique=True),
    "channel": pa.Column(str,
                         checks=pa.Check(lambda s: s.isin(criteria['allowed_channels'])),
                         required=True,
                         unique=True),
    "dye": pa.Column(str,
                     checks=pa.Check(lambda s: s.isin(criteria['allowed_dyes'])),
                     required=True,
                     unique=True),
    "ladder": pa.Column(bool,
                        checks=pa.Check(lambda b: sum(b) == criteria['ladder_channel_count']),
                        required=True,
                        unique=False),
    "desc": pa.Column(str,
                      checks=pa.Check.str_matches(criteria['regex_desc']),
                      required=False,
                      nullable=True,
                      unique=False)
})


def validator(indata):
    method_names = indata.keys()
    if len(method_names) != 1:
        raise ValueError

    df = pd.json_normalize(indata[method_names[0]])
    schema.validate(df)
    return True
