import matplotlib


plt_params_dict = {'font.sans-serif': 'Arial',
                   'mathtext.default': 'regular'}

def set_styles(plt, matplotlib):
    """
    Set fonts and default style
    """

    try:
        plt.style.use('default')        
        for param, value in plt_params_dict.items():
            print(param, value)
            plt.rcParams.update({param: value})
    except:
        print("""
            Before running set_styles(), you must:

            import matplotlib.pyplot as plt
            import matplotlib
            """)