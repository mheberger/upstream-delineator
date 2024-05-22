import psycopg2
from sqlalchemy import create_engine

host = "localhost"
database = "basins"
user = "postgres"
password = "dbpw"
port = 5432

# Connect to your PostgreSQL database
connection = psycopg2.connect(
    host=host,
    database=database,
    user=user,
    password=password
)
cursor = connection.cursor()
# Create a SQLAlchemy engine
engine = create_engine('postgresql://{}:{}@{}:{}/{}'.format(user, password, host, port, database))


def db_close():
    """
    Commit, close the cursor and connections
    """
    connection.commit()
    cursor.close()
    connection.close()
    engine.dispose()