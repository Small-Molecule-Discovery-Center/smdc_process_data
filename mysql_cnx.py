"""SQLAlchemy connections configured at call time from environment variables."""
import os
from sqlalchemy import create_engine
from sqlalchemy.engine import URL


def create_mysql_cnx(host="localhost", db=None):
    """Use mysql80 credentials locally and mysql credentials for remote hosts."""
    host = host or os.environ.get("hits_host", "localhost")
    db = db or os.environ.get("hits_db", "QB3_HiTS")
    local = host in {"localhost", "127.0.0.1", "::1"}
    password_key = "mysql80_pwd" if local else "mysql_pwd"
    port_key = "mysql80_port" if local else "mysql_port"
    url = URL.create(
        "mysql+pymysql", username=os.environ["mysql_user"],
        password=os.environ[password_key], host=host,
        port=int(os.environ[port_key]), database=db,
    )
    return create_engine(url)
