from pydantic import BaseModel, ValidationError, validator, root_validator
from typing import List, Literal, Optional
import re
import yaml
from yaml.loader import SafeLoader
import portion as P
import itertools
import logging

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
    """
    A Pydantic class for channel.
    ...

    Attributes
    ----------
    channel_name : str channel name
    channel : Literal of allowed channels
    dye : Literal of allowed dyes
    desc : Optional string
    alleles: List of Allele objects

    Methods
    -------
    channel_name_test(cls, value):
        Tests if value is matches the supplied regex

    unique_allele_names(cls, values):
        Validates that allele objects are unique

    non_overlapping_allele_intervals(cls, values):
        Validates that allele positions are non-overlapping

    """

    channel_name: str
    channel: Literal[tuple(c['allowed_channels'])]
    dye: Literal[tuple(c['allowed_dyes'])]
    ladder: bool
    desc: Optional[str]
    alleles: List[Allele]

    @validator('channel_name')
    def channel_name_test(cls, v):
        if not re.search(c['channel_name_regex'], v):
            raise ValueError('Channel name invalid')
        return v

    @root_validator()
    def unique_allele_names(cls, values):

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
    """
    A Pydantic class for method.
    ...

    Attributes
    ----------
    name : str
        method name
    channels : List of Channel objects
        channel

    Methods
    -------
    name_test(cls, value):
        Tests if value is matches the supplied regex
    """
    name: str
    channels: List[Channel]

    @validator('name')
    def name_test(cls, v):
        """
        Test if value is matches the supplied regex
        :param v: name field
        :return: name field if match is successful
        """

        if not re.search(c['method_name_regex'], v):
            raise ValueError('Method name invalid')
        return v

def create_method_obj(method: dict, logger: logging.Logger = None) -> Method:
    """

    :param method:
    :param logger:
    :return:
    """
    if logger:
        try:
            method_obj = Method(**method)
            logger.info(f'Method {method_obj.name} validated')
            return method_obj
        except Exception as e:
            logger.error(f"Validation of method failed")
            for row in str(e).split('\n'):
                logger.error(f"Validation error:     {row}")

            exit(1)

    else:
        method_obj = Method(**method)
        return method_obj

