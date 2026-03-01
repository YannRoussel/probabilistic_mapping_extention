# download e-features from EFel

## BBP

### collection of cells
aws s3 sync s3://openbluebrain/Model_Data/Brain_atlas/Mouse/ttypes_efeatures/releases/18122023/sscx/ ./sscx --no-sign-request

### cell_wise
aws s3 sync s3://openbluebrain/Model_Data/Brain_atlas/Mouse/ttypes_efeatures/releases/20122023/sscx_cellwise/ ./sscx_cellwise --no-sign-request


###cell_wise e-types
aws s3 sync s3://openbluebrain/Model_Data/Brain_atlas/Mouse/ttypes_efeatures/releases/20122023/sscx_cellwise_petilla_etypes.py ./sscx_cellwise_petilla_etypes.py --no-sign-request


## Scala

### 34 degrees 
aws s3 sync s3://openbluebrain/Model_Data/Brain_atlas/Mouse/ttypes_efeatures/releases/18122023/scala_34/ ./scala_34 --no-sign-request


### RT
aws s3 sync s3://openbluebrain/Model_Data/Brain_atlas/Mouse/ttypes_efeatures/releases/19122023/scala_roomtemp/ ./scala_roomtemp --no-sign-request


## Gouwens
aws s3 sync s3://openbluebrain/Model_Data/Brain_atlas/Mouse/ttypes_efeatures/releases/20122023/allen_34_VCtx_partial/ ./allen_34_VCtx_partial --no-sign-request
