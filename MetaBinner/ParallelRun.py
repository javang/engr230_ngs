
import multiprocessing

class ParallelRun:
    """ Class to run a function in parallel """
    def __init__(self):
        self.failed_processes = set()


    def run(self, func, argslist, kwargs, error_check=True):
        """
            Runs the function using all the arguemnts in argslist
            @param func A function.
            @param argslist A list of tuples. Each tuple contains the
              parameters for the function.
            @param kwargs. This list a dictionar with the values of the
              parameters for the function that are not variable.

            @return Returs a list with the values of the function
            @note The function get_failed_processes allows to recover
              the values of the parameters that failed
        """
        pool = multiprocessing.Pool(multiprocessing.cpu_count()-1)
        handles = []
        for args in zip(argslist):
            handle = pool.apply_async(func,args=args, kwds=kwargs)
            handles.append(handle)
        pool.close()
        pool.join()
        results = []
        for i, r in enumerate(handles):
            if not error_check:
                results.append(r.get())
            else:
                try:
                    results.append(r.get())
                except Exception as e:
                    self.failed_processes.add((argslist[i], e))
        return results

    def get_failed_processes(self):
        """ Returns the arguments used for the processes that failed

            @return A list of tuples of the form: (arguments, failure message)
        """
        return self.failed_processes


def hola(name, word="hola", sign="?"):
    return word + " " + name + sign


if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser(
        description="""Apply machine learning algorithms to the k-mer spectrums
                        of the scaffolds from the metagenome.
                    """)
    parser.add_argument("fn_database",
                    help="Datbabase formed by the files provided by the IMG/M for a metagenome. " \
                    "It must be created using the create_database.py script.")
    args = parser.parse_args()

    names = ["javi", "luis", "nacho"]
    mpr = ParallelRun()
    results = mpr.run(hola, names, {"word": "hi", "sign": " !"})
    for r in results:
        print r
    for p in mpr.get_failed_processes():
        print p