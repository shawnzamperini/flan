# distutils: language=c++

# Class to help with reading openadas files.
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import griddata


class OpenADAS:

    def __init__(self):
        pass

    @staticmethod
    def read_rate_coef_unres(path):
        """
        Read in the rate coefficients from a given unresolved file. This could be an ACD, SCD, XCD, QCD or CCD file. A
        description of the file formatting convention is given here: https://open.adas.ac.uk/man/appxa-11.pdf.

        path (str): Path to the file.
        """

        with open(path, "r") as f:

            # Get header info.
            header = f.readline().split()
            atomic_number = int(header[0])
            ndens = int(header[1])
            ntemp = int(header[2])
            charge_low = int(header[3])
            charge_high = int(header[4])
            element = header[5][1:].title()
            # year = int(header[-1])
            ncharges = charge_high - charge_low + 1
            print("ADAS: Reading in data for {} Z={}...".format(element, atomic_number))
            f.readline()

            # Setup arrays.
            ne = np.zeros(ndens)
            te = np.zeros(ntemp)

            # Read in densities (actually log10(ne) here).
            count = 0
            while count < ndens:
                line = f.readline().split()
                for i in range(0, len(line)):
                    ne[count + i] = line[i]
                count += len(line)
                if count == ndens:
                    break
                elif count > ndens:
                    print("Error! count > ndens ({} > {})".format(count, ndens))
                    break

            # Again for temperatures.
            count = 0
            while count < ntemp:
                line = f.readline().split()
                for i in range(0, len(line)):
                    te[count + i] = line[i]
                count += len(line)
                if count == ntemp:
                    break
                elif count > ntemp:
                    print("Error! count > ntemp ({} > {})".format(count, ntemp))
                    break

            # Assemble a DataFrame to be filled out with the rate coefficients.
            multi_index = pd.MultiIndex.from_product([list(range(1, ncharges+1)), ne], names=["charge", "density"])
            df = pd.DataFrame(index=multi_index, columns=te, dtype=float)

            # Go through one charge state at a time.
            for charge in range(1, ncharges+1):
                f.readline()

                # The file prints out the rate coefficients in groups of same Te's.
                for t in te:
                    count = 0
                    while count < ndens:
                        line = f.readline().split()
                        for val in line:
                            df.loc[charge, ne[count]][t] = float(val)
                            count += 1
                        if count == ntemp:
                            break
                        elif count > ntemp:
                            print("Error! count > ntemp ({} > {}) for charge {}".format(count, ntemp, charge))
                            break

            return df

    @staticmethod
    def get_rate_coef(df, te, ne, charge=1):
        """
        Using the DataFrame returned from read_rate_coef, return the rate coefficient at a given ne, Te by 2D linearly
        interpolating between the (log10) ne, Te values.

        Input
        te: Te to get rate coefficient at in eV.
        ne: ne to get rate coefficient at in m-3.

        Output
        Returns the rate coefficient in units of m3/s
        """

        tes = df.columns.values
        nes = df.loc[charge].index.values
        rates = df.loc[charge].values

        # griddata is the fastest 2D interpolater that I know of.
        X, Y = np.meshgrid(tes, nes)
        points = np.array((X.flatten(), Y.flatten())).T
        xi = (np.log10(te), np.log10(ne * 1e-6))  # ne in the data is cm-3, hence the 1e-6.
        return np.power(10, griddata(points, rates.flatten(), xi)) * 1e-6

    @staticmethod
    def plot_rates(df, charge=1):
        """
        Whatever rates were loaded in from read_rate_coef, plot them here. Only tested so far with unresolved rate
        coefficients.

        df: DataFrame returned from read_rate_coef_unres.
        """

        tes = df.columns.values
        nes = df.loc[charge].index.values

        fig, ax1 = plt.subplots(figsize=(10, 8))

        for ne in nes:
            rates = df.loc[charge, ne].values
            ax1.plot(np.power(10, tes), np.power(10, rates), label="{:.2e}".format(np.power(10, ne)))

        ax1.grid()
        ax1.legend(ncol=2)
        ax1.set_xscale("log")
        ax1.set_yscale("log")
        ax1.set_xlabel("Te (eV)")
        ax1.set_ylabel(r"Rate Coefficient $\mathdefault{cm^3/s}$")
        fig.tight_layout()
        fig.show()

