from logging import getLogger
from os.path import splitext, basename

import colored_logging as cl
from .ECOSTRESS_granule import ECOSTRESSGranule
from .L1_RAD import L1CGRAD, L1CTRAD
from .L1B_GEO import L1BGEO
from .L2_LSTE import L2GLSTE, L2LSTE, L2TLSTE, ECOv002L2LSTE
from .L2_STARS import L2TSTARS
from .L3_ET_ALEXI import L3ETALEXI, L3TETALEXI, L3GETALEXI
from .L3_ET_PTJPL import L3TETPTJPL, L3TMET, L3TSM, L3TSEB, L3GETPTJPL, L3GMETPTJPL, L3GSMPTJPL, L3GSEBPTJPL, L3TJET, L3ETPTJPL
from .L3_JET import L3GJET
from .L4_ESI import L4ESI, L4TESI, L4GESI
from .L4_ESI_ALEXI import L4ESIALEXI, L4TESIALEXI, L4GESIALEXI
from .L4_WUE import ECOv001L4WUEPTJPL, L4TWUE, L4GWUE

cl.configure()
logger = getLogger(__name__)


def open_granule(
        filename: str,
        L2_CLOUD_filename: str = None,
        L1B_GEO_filename: str = None,
        **kwargs) -> ECOSTRESSGranule:
    filename_base = splitext(basename(filename))[0]

    if filename_base.startswith("ECOv002_L1CG_RAD"):
        logger.info(f"loading Collection 2 L2G RAD: {cl.file(filename)}")
        return L1CGRAD(filename, **kwargs)
    elif filename_base.startswith("ECOv002_L1CT_RAD"):
        logger.info(f"loading Collection 2 L2T RAD: {cl.file(filename)}")
        return L1CTRAD(filename, **kwargs)
    elif filename_base.startswith("ECOv002_L2G_LSTE"):
        logger.info(f"loading Collection 2 L2G LSTE: {cl.file(filename)}")
        return L2GLSTE(filename, **kwargs)
    elif filename_base.startswith("ECOv002_L2T_LSTE"):
        logger.info(f"loading Collection 2 L2T LSTE: {cl.file(filename)}")
        return L2TLSTE(filename, **kwargs)
    elif filename_base.startswith("ECOv002_L2T_STARS"):
        logger.info(f"loading Collection 2 L2T STARS: {cl.file(filename)}")
        return L2TSTARS(filename, **kwargs)
    elif filename_base.startswith("ECOv002_L3T_JET"):
        logger.info(f"loading Collection 2 L3T JET: {cl.file(filename)}")
        return L3TJET(filename, **kwargs)
    elif filename_base.startswith("ECOv002_L3G_JET"):
        logger.info(f"loading Collection 2 L3G JET: {cl.file(filename)}")
        return L3GJET(filename, **kwargs)
    elif filename_base.startswith("ECOv002_L3T_ET_ALEXI"):
        logger.info(f"loading Collection 2 L3T ET ALEXI: {cl.file(filename)}")
        return L3TETALEXI(filename, **kwargs)
    elif filename_base.startswith("ECOv002_L3G_ET_ALEXI"):
        logger.info(f"loading Collection 2 L3G ET ALEXI: {cl.file(filename)}")
        return L3GETALEXI(filename, **kwargs)
    elif filename_base.startswith("ECOv002_L3T_MET"):
        logger.info(f"loading Collection 2 L3T MET: {cl.file(filename)}")
        return L3TMET(filename, **kwargs)
    elif filename_base.startswith("ECOv002_L3G_MET"):
        logger.info(f"loading Collection 2 L3G MET: {cl.file(filename)}")
        return L3GMETPTJPL(filename, **kwargs)
    elif filename_base.startswith("ECOv002_L3T_SM"):
        logger.info(f"loading Collection 2 L3T SM: {cl.file(filename)}")
        return L3TSM(filename, **kwargs)
    elif filename_base.startswith("ECOv002_L3G_SM"):
        logger.info(f"loading Collection 2 L3G SM: {cl.file(filename)}")
        return L3GSMPTJPL(filename, **kwargs)
    elif filename_base.startswith("ECOv002_L3T_SEB"):
        logger.info(f"loading Collection 2 L3T SEB: {cl.file(filename)}")
        return L3TSEB(filename, **kwargs)
    elif filename_base.startswith("ECOv002_L3G_SEB"):
        logger.info(f"loading Collection 2 L3G SEB: {cl.file(filename)}")
        return L3GSEBPTJPL(filename, **kwargs)
    elif filename_base.startswith("ECOv002_L4T_ESI_ALEXI"):
        logger.info(f"loading Collection 2 L4T ESI ALEXI: {cl.file(filename)}")
        return L4TESIALEXI(filename, **kwargs)
    elif filename_base.startswith("ECOv002_L4G_ESI_ALEXI"):
        logger.info(f"loading Collection 2 L4G ESI ALEXI: {cl.file(filename)}")
        return L4GESIALEXI(filename, **kwargs)
    elif filename_base.startswith("ECOv002_L4T_ESI"):
        logger.info(f"loading Collection 2 L4T ESI: {cl.file(filename)}")
        return L4TESI(filename, **kwargs)
    elif filename_base.startswith("ECOv002_L4G_ESI"):
        logger.info(f"loading Collection 2 L4G ESI: {cl.file(filename)}")
        return L4GESI(filename, **kwargs)
    elif filename_base.startswith("ECOv002_L4T_WUE"):
        logger.info(f"loading Collection 2 L4T WUE: {cl.file(filename)}")
        return L4TWUE(filename, **kwargs)
    elif filename_base.startswith("ECOv002_L4G_WUE"):
        logger.info(f"loading Collection 2 L4G WUE: {cl.file(filename)}")
        return L4GWUE(filename, **kwargs)
    elif filename_base.startswith("ECOSTRESS_L1B_GEO"):
        logger.info(f"loading Collection 1 L1B GEO: {cl.file(filename)}")
        return L1BGEO(filename, **kwargs)
    elif filename_base.startswith("ECOSTRESS_L2_LSTE"):
        logger.info(f"loading Collection 1 L2 LSTE: {cl.file(filename)}")
    
        return L2LSTE(
            L2_LSTE_filename=filename, 
            L2_CLOUD_filename=L2_CLOUD_filename,
            L1B_GEO_filename=L1B_GEO_filename,
            **kwargs
        )
    
    elif filename_base.startswith("ECOSTRESS_L3_ET_PT-JPL"):
        logger.info(f"loading Collection 1 L3 ET PT-JPL: {cl.file(filename)}")
        
        return L3ETPTJPL(
            L3_ET_PTJPL_filename=filename, 
            L2_CLOUD_filename=L2_CLOUD_filename,
            L1B_GEO_filename=L1B_GEO_filename,
            **kwargs
        )
    
    elif filename_base.startswith("ECOSTRESS_L3_ET_ALEXI"):
        logger.info(f"loading Collection 1 L3 ET ALEXI: {cl.file(filename)}")
        
        return L3ETALEXI(
            L3_ET_ALEXI_filename=filename, 
            L2_CLOUD_filename=L2_CLOUD_filename,
            L1B_GEO_filename=L1B_GEO_filename,
            **kwargs
        )
    
    elif filename_base.startswith("ECOSTRESS_L4_ESI_PT-JPL"):
        logger.info(f"loading Collection 1 L4 ESI PT-JPL: {cl.file(filename)}")
        
        return L4ESI(
            L4_ESI_filename=filename, 
            L2_CLOUD_filename=L2_CLOUD_filename,
            L1B_GEO_filename=L1B_GEO_filename,
            **kwargs
        )
    
    elif filename_base.startswith("ECOSTRESS_L4_ESI_ALEXI"):
        logger.info(f"loading Collection 1 L4 ESI ALEXI: {cl.file(filename)}")
        
        return L4ESIALEXI(
            L4_ESI_ALEXI_filename=filename, 
            L2_CLOUD_filename=L2_CLOUD_filename,
            L1B_GEO_filename=L1B_GEO_filename,
            **kwargs
        )
    
    elif filename_base.startswith("ECOSTRESS_L4_WUE"):
        logger.info(f"loading Collection 1 L4 WUE: {cl.file(filename)}")

        return ECOv001L4WUEPTJPL(
            L4_WUE_filename=filename,
            L2_CLOUD_filename=L2_CLOUD_filename,
            L1B_GEO_filename=L1B_GEO_filename,
            **kwargs
        )
    else:
        logger.warning(f"unrecognized file: {filename}")
