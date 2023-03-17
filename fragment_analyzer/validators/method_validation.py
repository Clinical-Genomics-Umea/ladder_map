from pydantic import BaseModel, ValidationError, validator, root_validator
from typing import List, Literal
import re
import yaml
from yaml.loader import SafeLoader
import portion as P
import itertools

with open('fragment_analyzer/config/method_validation_criteria.yaml') as f:
    c = yaml.load(f, Loader=SafeLoader)

#
# def root_unique_validator(field):
#     def validator(cls, values):
#         root_values = values.get("__root__")
#         value_set = set()
#         for value in root_values:
#             if value[field] in value_set:
#                 raise ValueError(f"Duplicate {field}")
#             else:
#                 value_set.add(value[field])
#         return values
#
#     return root_validator(pre=True, allow_reuse=True)(validator)


class Allele(BaseModel):
    allele_name: str
    start: int
    end: int

    @validator('allele_name')
    def allele_test(cls, v):
        if not re.search(c['allele_name_regex'], v):
            raise ValueError('Allele name invalid')
        return v.title()


class Channel(BaseModel):
    channel_name: str
    channel: Literal[tuple(c['allowed_channels'])]
    dye: Literal[tuple(c['allowed_dyes'])]
    ladder: bool
    desc: str
    alleles: List[Allele]

    @validator('channel_name')
    def name_test(cls, v):
        if not re.search(c['channel_name_regex'], v):
            raise ValueError('Channel name invalid')
        return v.title()

    # @validator('channel')
    # def name_test(cls, v):
    #     if v not in c['allowed_channels']:
    #         raise ValueError('Channel invalid')
    #     return v.title()

    # @validator('dye')
    # def name_test(cls, v):
    #     if v not in c['allowed_dyes']:
    #         raise ValueError('Dye invalid')
    #     return v.title()

    # @validator('desc')
    # def name_test(cls, v):
    #     if not re.search(c['desc_regex'], v):
    #         raise ValueError('Description invalid')
    #     return v.title()

    @root_validator()
    def unique_allele_values(cls, values):

        uq_fields = c['alleles_unique_fields']
        alleles_list = values.get('alleles')

        for uf in uq_fields:
            a_vals = [getattr(a, uf) for a in alleles_list]

            a_vals_set = set(a_vals)
            if len(a_vals) != len(a_vals_set):
                 raise ValueError(f'Duplicates in allele field {uf}')

        return values

    @root_validator()
    def non_overlapping_allele_intervals(cls, values):

        alleles_list = values.get('alleles')
        intervals = [P.closed(a.start, a.end) for a in alleles_list]

        ipairs = itertools.combinations(intervals, 2)

        for i1, i2 in ipairs:
            res = i1 & i2
            if not res == P.empty():
                raise ValueError(f'Alleles {i1} and {i2} overlap')

        return values


class Method(BaseModel):
    name: str
    channels: List[Channel]

    @validator('name')
    def name_test(cls, v):
        if not re.search(c['method_name_regex'], v):
            raise ValueError('Method name invalid')
        return v.title()
    #
    # @root_validator(pre=True)
    # def unique_values(cls, values):
    #     root_values = values.get('channels')
    #
    #     fields = c['channels_unique_fields']
    #     value_sets = {}
    #     for f in fields:
    #         value_sets[f] = set()
    #
    #     for value in root_values:
    #         for f in fields:
    #             if value[f] in value_sets[f]:
    #                 raise ValueError(f'Duplicate in field {f}')
    #             else:
    #                 value_sets[f].add(value[f])
    #
    #     return values

    # @root_validator()
    # def unique_channel_values(cls, values):
    #
    #     uq_fields = c['channels_unique_fields']
    #     channel_list = values.get('channels')
    #
    #     for uf in uq_fields:
    #         ch_vals = [getattr(ch, uf) for ch in channel_list]
    #
    #         ch_vals_set = set(ch_vals)
    #         if len(ch_vals) != len(ch_vals_set):
    #             raise ValueError(f'Duplicates in channel field {uf}')
    #
    #     return values

def validate_method(method: dict):
    return Method(**method)


