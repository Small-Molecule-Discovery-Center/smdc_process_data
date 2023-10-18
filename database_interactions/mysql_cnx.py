# import mysql.connector
from sqlalchemy import create_engine
import os


###############################################################################
def create_mysql_cnx(host=os.environ['hits_host'], db=os.environ['hits_db']):
    """ A helper function to create a sql alchemy engine by pulling connection 
    info from environment.
    """
    cnx = create_engine(f"mysql+pymysql://{os.environ['mysql_user']}:{os.environ['mysql_pwd']}@{host}:{os.environ['mysql_port']}/{db}")
    return cnx