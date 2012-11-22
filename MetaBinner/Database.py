
import sqlite3 as sqlite
import os
import csv
import logging

log = logging.getLogger("Database")

class Database3(sqlite.Connection):
    """ Class to manage a SQL database built with sqlite3 """

    name_to_type = {'INT':int, 'DOUBLE':float, 'VARCHAR(10)':str}
    type_to_name = {int:'INT', float:'DOUBLE', str:'VARCHAR(10)'}

    def __init__(self, fn_database, overwrite=False):
        """ Creates a database by simply connecting to the file

            @param overwrite If true, the database is deleted and created again
        """
        if overwrite and os.path.exists(fn_database):
            log.info("Creating database")
            os.remove(fn_database)
        sqlite.Connection.__init__(self, fn_database)
        self.row_factory = sqlite.Row

    def drop_table(self, table_name):
        """ Delete a table if it exists """
        self.execute("DROP TABLE IF EXISTS %s" % (table_name))
        self.commit()

    def retrieve_data(self,sql_command):
        """ Retrieves data from the database using the sql_command
        returns the records as a list of tuples"""
        return self.execute(sql_command).fetchall()

    def get_tables_names(self):
        """ Return the name of all the tables in the database """
        sql_command = """ SELECT tbl_name FROM sqlite_master """
        data = self.retrieve_data(sql_command)
        return frozenset([d[0] for d in data])

    def create_table(self, table_name, column_names, column_types):
        """ Creates a table.
            @param table_name The name of the table_name
            @param column_names List with the names of all the columns of the
             table
            @param column_type Type of all the values in the database. They
            are python types.
            Example: create_database("table",["number","string"],[int,str])
        """
        log.info("Creating table %s",table_name)
        sql_command = "CREATE TABLE %s (" % (table_name)
        for name, data_type in zip(column_names, column_types):
            sql_command += "{0} {1},".format(name, self.type_to_name[data_type])
        # replace last comma for a parenthesis
        self.execute(sql_command[0:-1] + ")")
        self.commit()

    def store_data(self,table_name,data):
        """ Inserts information in a given table of the database.
        The info must be a list of tuples containing as many values
        as columns in the table
            Conversion of values is done AUTOMATICALLY after checking the
            types stored in the table
        """
        if len(data) == 0:
            log.warning("Inserting empty data")
            return
        n = len(data[0]) # number of columns for each row inserted
        tuple_format="("+"?,"*(n-1)+"?)"
        sql_command="INSERT INTO %s VALUES %s " % (table_name, tuple_format)
        # Fill the table with the info in the tuples
        types = self.get_table_types(table_name)
#        log.debug("Storing types: %s", types)
        for i in xrange(len(data)):
            data[i] = [apply_type(d) for d,apply_type in zip(data[i], types)]
        self.executemany(sql_command, data)
        self.commit()

    def get_table_types(self, name):
        """ Retuns all the types in a table
        """
        sql_command = "PRAGMA table_info(%s)" % name
        info = self.execute(sql_command).fetchall()
        return [self.name_to_type[row[2]] for row in info]


    def get_table_column_names(self, name):
        """ Get the names of the columns for a given table
            @param name The name of the  table
        """
        cursor = self.execute("PRAGMA table_info(%s)" % name)
        return [ row[1] for row in cursor.fetchall()]

    def add_column(self,table, column, data_type):
        """ Add a column to a table
            @param column - the name of the column.
            @param data_type - the type: int, float, str
        """
        sql_command = "ALTER TABLE %s ADD %s %s" % (table, column, self.type_to_name[data_type])
        self.execute(sql_command)

    def add_columns(self, table, names, types, check=True):
        """ Add columns to the database.

            @param check If check=True, columns with names
            already in the database are skipped. If check=False no check
            is done and trying to add a column that already exists will
            raise and exception
        """
        if check:
            col_names = self.get_table_column_names(table)
            for name, dtype in zip(names, types):
                if name not in col_names:
                    self.add_column(table, name, dtype)
        else:
            for name, dtype in zip(names, types):
                self.add_column(table, name, dtype)



    def drop_view(self,view_name):
        """ Removes a view from the database """
        self.execute('DROP VIEW %s' % view_name)


    def drop_columns(self, table, columns):

        cnames = self.get_table_column_names(table)
        map(cnames.remove,columns)
        names_txt = ", ".join(cnames)
        sql_command = [
        "CREATE TEMPORARY TABLE backup(%s);" % names_txt,
        "INSERT INTO backup SELECT %s FROM %s" % (names_txt, table),
        "DROP TABLE %s;" % table,
        "CREATE TABLE %s(%s);" % (table, names_txt),
        "INSERT INTO %s SELECT * FROM backup;" % table,
        "DROP TABLE backup;",
        ]
        map(self.execute,sql_command)


    def update_data(self, table_name,
                    updated_fields,
                    updated_values,
                    condition_fields,
                    condition_values):
        """ updates the register in the table identified by the condition
            values for the condition fields
        """
        s = ["{0}={1}".format(f, v) for f,v in zip(updated_fields, updated_values)]
        update_str = ",".join(s)
        s = ["{0}={1}".format(f, v) for f,v in zip(condition_fields, condition_values)]
        condition_str = " AND ".join(s)
        sql_command = """ UPDATE {0}
                          SET {1}
                          WHERE {2}
                      """.format(table_name, update_str, condition_str)
        log.debug("Updating %s: %s",table_name, sql_command)
        self.execute(sql_command)
        self.commit()




    def create_view(self,view_name,table_name,
                                condition_fields, condition_values):
        """ creates a view of the given table where the values are selected
            using the condition values. See the help for update_data()
        """
        try: # if this fails is because the view does not exist
           self.drop_view(view_name)
        except:
            pass
        s = ["{0}={1}".format(f, v) for f,v in zip(condition_fields, condition_values)]
        condition_str = " AND ".join(s)
        sql_command = """ CREATE VIEW {0}
                          AS SELECT * FROM {1}
                          WHERE {2} """.format(view_name, table_name, condition_str)
        log.info("Creating view %s", sql_command)
        self.execute(sql_command)


    def show(self):
        """ Show the names of the tables and columns in the database in a friendly format
        """
        for t in self.get_tables_names():
            print t
            for c in self.get_table_column_names(t):
                print "   |_ {0}".format(c)

class Database2:


    def create_view_of_best_records(self, view_name, table_name, orderby, n_records):
        try: # if this fails is because the view already exist
           self.drop_view(view_name)
        except:
            pass
        sql_command = """CREATE VIEW %s AS SELECT * FROM %s
                         ORDER BY %s ASC LIMIT %d """ % (view_name, table_name, orderby, n_records)
        log.info("Creating view %s", sql_command)
        self.cursor.execute(sql_command)


    def get_table(self, table_name, fields=False, orderby=False):
        """ Returns th fields requested from the table """
        fields = self.get_fields_string(fields)
        sql_command = "SELECT %s  FROM %s " % (fields, table_name)
        if orderby:
            sql_command += " ORDER BY %s ASC" % orderby
        data = self.retrieve_data(sql_command)
        return data

    def get_fields_string(self, fields, field_delim=","):
        if fields:
            return field_delim.join(fields)
        return "*"



    def select_table(self):
        """
            Prompt for tables so the user can choose one
        """
        table_name = ""
        self.check_if_is_connected()
        tables = self.get_tables_names()
        for t in tables:
            say = ''
            while say not in ('n','y'):
                say = raw_input("Use table %s (y/n) " % t)
            if say == 'y':
                table_name = t
                columns = self.get_table_column_names(t)
                break
        return table_name, columns




def print_data(data, delimiter=" "):
    """ Prints the data recovered from a database """
    for row in data:
        line = delimiter.join([str(x) for x in row])
        print line

def write_data(data,output_file,delimiter=" "):
    """writes data to a file. The output file is expected to be a python
    file object """
    w = csv.writer(output_file, delimiter=delimiter)
    map(w.writerow, data)


def read_data(fn_database, sql_command):
    db = Database3(fn_database)
    data = db.retrieve_data(sql_command)
    db.close()
    return data

def get_sorting_indices(l):
    """ Return indices that sort the list l """
    pairs = [(element,i) for i,element in enumerate(l)]
    pairs.sort()
    indices = [p[1] for p in pairs]
    return indices

def merge_databases(fns, fn_output, tbl):
    """
       Reads a table from a set of database files into a single file
       Makes sure to reorder all column names if neccesary before merging
    """
    # Get names and types of the columns from first database file
    db = Database3(fns[0])
    names = db.get_table_column_names(tbl)
    types = db.get_table_types(tbl)
    indices = get_sorting_indices(names)
    sorted_names = [ names[i] for i in indices]
    sorted_types = [ types[i] for i in indices]
    log.info("Merging databases. Saving to %s", fn_output)
    out_db = Database3(fn_output, overwrite=True)
    out_db.create_table(tbl, sorted_names, sorted_types)
    for fn in fns:
        log.debug("Reading %s",fn)
        names = db.get_table_column_names(tbl)
        names.sort()
        they_are_sorted = ",".join(names)
        log.debug("Retrieving %s", they_are_sorted)
        sql_command = "SELECT %s FROM %s" % (they_are_sorted, tbl)
        data = db.retrieve_data(sql_command)
        out_db.store_data(tbl, data)
        db.close()
    out_db.close()

