"""A set of functions to generally help create new entries in the metadata sections of HiTS DB."""

from datetime import datetime
from sqlalchemy import Table, MetaData, select, text, insert, func
import pandas as pd
import numpy as np

def string_now():
    """Create a date / time formatted as a string for MySQL.
    Returns:
        A string version of the date and time.
    """
    return datetime.now().strftime('%Y-%m-%d %H:%M:%S')


def check_unique(new_df, table, unique_col, engine):
    """Check a column in a HiTS table for uniqueness.
    Args:
        new_df (pd.DataFrame): a dataframe with columns named for the database table and rows containing new data to be inserted.
        table (sqlalchemy.Table): a table from the database to check against.
        unique_col (str): the string name of the column to check. Must be in both new_df and table.
        engine (sqlalchemy connection): the engine you want to use to connect to the database.
    Returns:
        A boolean; True if the entries are unique in the column and false otherwise.
    """
    already_exists=False
    stmt=select(table).where(table.c[unique_col].in_(tuple(new_df[unique_col].tolist())))
    with engine.connect() as conn:
        result=conn.execute(stmt)
        already_exists=result.rowcount>0
    
    if already_exists:
        print([x.name for x in table.columns])
        for i, row in enumerate(result):
            print(row)
            # assert i<0, 
            print(f"\nWarning: this {unique_col} already exists")
            return False
    else:
        print(f'OK - these {unique_col} values are unique.')
        return True


def count_rows(table, engine):
    """Count rows in a table to gauge size before using with Pandas.
    Args:
        table (sqlalchemy Table): the table to count.
        engine (sqlalchemy connection): the engine to use to connect to the database.
    Returns:
        Number of rows based on the primary key of the table.
    """
    print(f"Primary key: {table.primary_key.c[0].name}")
    with engine.connect() as conn:
        stmt=select(func.count(table.c[table.primary_key.c[0].name]))
        result=conn.execute(stmt)
    return result.all()