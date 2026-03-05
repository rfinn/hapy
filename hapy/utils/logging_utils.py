from pathlib import Path
import logging


def setup_logging(outdir=None, tag=None, script_name=None, level=logging.INFO):
    """
    Configure console + file logging for pipeline scripts.

    Parameters
    ----------
    outdir : str or Path
        Base output directory where logs/ will be created.
    tag : str
        Identifier for this run (e.g. image_id or galaxy tag).
    script_name : str
        Name of the script (e.g. 'cutouts', 'analysis').
    level : int
        Logging level (default INFO)

    Returns
    -------
    logger : logging.Logger
    """

    # configure root logger once
    logging.basicConfig(
        level=level,
        format="%(asctime)s %(levelname)s [pid=%(process)d] %(name)s: %(message)s",
    )

    root = logging.getLogger()

    # build logfile path
    if outdir:
        outdir = Path(outdir)
        logdir = outdir / "logs"
        logdir.mkdir(parents=True, exist_ok=True)

        name_parts = []
        if tag:
            name_parts.append(tag)
        if script_name:
            name_parts.append(script_name)

        logfile = logdir / (".".join(name_parts) + ".log")

        # avoid duplicate handlers
        if not any(isinstance(h, logging.FileHandler) for h in root.handlers):
            fh = logging.FileHandler(logfile)
            fh.setLevel(level)
            fh.setFormatter(logging.Formatter(
                "%(asctime)s %(levelname)s [pid=%(process)d] %(name)s: %(message)s"
            ))
            root.addHandler(fh)

    # return a tagged logger
    if tag:
        logger_name = f"{script_name}.{tag}" if script_name else tag
    else:
        logger_name = script_name or __name__

    return logging.getLogger(logger_name)
