import random

import matplotlib.pyplot as plt
import seaborn


def compute_precision_recall_point(retrieved, relevant):
    '''
    Computes precision/recall for a set of items.

    Precision = | relevant intersect retrieved | / | retrieved|
    Recall = | relevant intersect retrieved | / | relevant |

    :param retrieved: Set of retrieved items
    :param relevant: Set of relevant items
    '''
    precision = float(len(
        relevant.intersection(retrieved))) / len(retrieved)

    recall = float(len(
        relevant.intersection(retrieved))) / len(relevant)

    return (precision, recall)


def compute_precision_recall_curve(ranked_items, relevant):
    '''
    Computes precision/recall for a list of ranked items.

    Precision = | relevant intersect retrieved | / | retrieved|
    Recall = | relevant intersect retrieved | / | relevant |

    :param ranked_items: List of sets of ranked items
    :param relevant: Set of relevant items

    '''

    retrieved = set()
    points = []
    points.append((1, 0))

    for items_set in ranked_items:
        retrieved = retrieved.union(items_set)

        precision = float(len(
            relevant.intersection(retrieved))) / len(retrieved)

        recall = float(len(
            relevant.intersection(retrieved))) / len(relevant)

        points.append(compute_precision_recall_point(retrieved, relevant))
        
    return points


def compute_precision_recall_curve_subsample_negatives(
        ranked_items, relevant, negatives, ratio=50):
    """
    negatives: set of all possible negatives
    ratio: subsample x * sizeof(negatives)
    """

    retrieved = set()
    points = []
    points.append((1, 0))

    subsample = None

    # Subsample negatives
    if len(negatives) <= len(relevant) * ratio:
        print("Warning: not enough negatives to meet desired ratio!")
        subsample = negatives
    else:
        neg_list = list(negatives)
        subsample = set(random.sample(negatives, len(relevant) * ratio))

    # Define universe of interactions we care about
    universe = relevant.union(subsample)

    # Run calculation, only caring about items in universe above 
    for items_set in ranked_items:
        for item in items_set:
            if item in universe:
                retrieved.add(item)
       
        if len(retrieved) > 0:
            precision = float(len(
                relevant.intersection(retrieved))) / len(retrieved)

            recall = float(len(
                relevant.intersection(retrieved))) / len(relevant)

            points.append(compute_precision_recall_point(retrieved, relevant))
        
    return points


def init_precision_recall_figure(): 
    fig, ax = plt.subplots()  
    ax.autoscale(False)
    ax.set_xlabel("Recall")
    ax.set_ylabel("Precision")
    ax.set_title("Precision Recall")
    ax.set_xbound(0, 1)
    ax.set_ybound(0, 1)

    return fig, ax


def plot_precision_recall_curve(ranked_items, relevant, ax=None):
    curve_points = compute_precision_recall_curve(ranked_items, relevant)

    x = [p[1] for p in curve_points]
    y = [p[0] for p in curve_points]

    if ax == None:
        fig, ax = init_precision_recall_figure()
        ax.plot(x, y)
        return fig, ax

    else:
        ax.plot(x, y)


def plot_precision_recall_curve_subsample_negatives(
        ranked_items, relevant, negatives, ratio=50, ax=None):
    curve_points = compute_precision_recall_curve_subsample_negatives(  
        ranked_items, relevant, negatives, ratio)

    x = [p[1] for p in curve_points]
    y = [p[0] for p in curve_points]

    if ax == None:
        fig, ax = init_precision_recall_figure()
        ax.plot(x, y)
        return fig, ax

    else:
        ax.plot(x, y)


if __name__=='__main__':
    main(sys.argv)
