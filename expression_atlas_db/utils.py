import numpy as np
import pandas as pd
from typing import List, Tuple, Dict
from pathlib import Path
import urllib
import xml.etree.ElementTree as ET
from io import StringIO
import time
import http
import boto3
import anndata as ad

from expression_atlas_db import settings

class GTFParser:
    """
    """
    genome_build = None
    genome_accession = None

    def __init__(self, gtf_path:str, drop_duplicate_id:bool=False):
        self.gtf_path = Path(gtf_path)
        self.gtf_name = self.gtf_path.stem
        self.parse_gtf(self.gtf_path, drop_duplicate_id=drop_duplicate_id)
    
    def parse_gtf(self, gtf_path:Path, drop_duplicate_id:bool=False) -> None:
        """Read GTF from file and create table structures.
        
        Args:
            gtf_path (str) path to gtf file. 
            drop_duplicate_id (bool) veliadb has duplicated gene/transcript ids. Drop and keep first.
        """
        transcripts = []
        genes = []
        
        with open(gtf_path, 'r') as f_in:
            transcript_count = 0
            gene_count = 0
            for line in f_in:
                if line.startswith('#'):
                    if line.startswith('#!genome-build '):
                            self.genome_build = line.split(' ')[-1].strip()
                    elif line.startswith('#!genome-build-accession '):
                        self.genome_accession = line.split(' ')[-1].strip()
                    continue
            
                fields, infos = self.parse_gtf_line(line.strip())

                if fields['label'] in ('gene','transcript',):
                    fields.update(infos)
                    if fields['label'] == 'transcript':
                        transcripts.append(fields)
                        transcript_count += 1

                    elif fields['label'] == 'gene':
                        genes.append(fields)
                        gene_count += 1

        if not self.genome_build:
            self.genome_build = gtf_path.stem.split('.gtf')[0]
                    
        transcript_df = pd.DataFrame(
                                transcripts,
                                columns=[
                                    'transcript_id',
                                    'gene_id',
                                    'veliadb_id',
                                    'gene_biotype',
                                    ],
                                )
        gene_df = pd.DataFrame(
                            genes,
                            columns=[
                                'gene_id',
                                'veliadb_id',
                                'gene_biotype',
                                ],
                            )
        
        transcript_df['veliadb_id'] = transcript_df['veliadb_id'].astype(int)
        gene_df['veliadb_id'] = gene_df['veliadb_id'].astype(int)
        transcript_df['assembly_id'] = self.genome_build
        gene_df['assembly_id'] = self.genome_build
                        
        if drop_duplicate_id:
            transcript_df.drop_duplicates('transcript_id', inplace=True)
            gene_df.drop_duplicates('gene_id', inplace=True)

        # if (transcript_df['transcript_id'].value_counts()[0] > 1) | \
        #     (transcript_df['veliadb_id'].value_counts()[0] > 1) | \
        #     (gene_df['gene_id'].value_counts()[0] > 1) | \
        #     (gene_df['veliadb_id'].value_counts()[0] > 1):
        #     raise ValueError('Duplicated gene/transcript/veliadb ids detected in gtf.')

        self.transcript_df = transcript_df.copy()
        self.gene_df = gene_df.copy()
              
    @staticmethod
    def parse_gtf_line(gtf_line: str, no_chr:bool=True) -> Tuple[Dict, Dict[str,str]]:
        """Parse gtf line and return fields.

        Args:
            gtf_line (str) raw gtf line
            no_chr (bool) flag to remove chr from chromosome names.

        Returns:
            (Tuple[Dict, Dict[str,str]]) tuple of dict of fields, dict of accessory line info.
        """
        chrom, source, label, start, end, score, strand, frame, info = gtf_line.split('\t')
        
        if no_chr:
            if chrom.startswith('chr'):
                chrom = chrom[3:]

            if chrom == 'M':
                chrom = 'MT'

        fields = {
                    'Chromosome': chrom,
                    'source': source, 
                    'label': label,
                    'Start': int(start),
                    'End': int(end),
                    'Score': score,
                    'Strand': strand,
                    'Frame': frame,
            }
        infos = {
                i.strip().split(' ')[0].strip(): i.strip().split(' ')[1].strip().replace('"','') \
                for i in info.strip(';\n').split(';')
            } 
        return fields, infos
    
class ExperimentParser:
    """
    """

    _h5_params = {}
    _s3_enabled = False
    _client = None

    def __init__(
                self, 
                velia_study_id:str, 
                exp_loc:Path=Path(settings.test_experiment_loc),
                ):
        self.velia_study_id = velia_study_id
        self.velia_study_loc = exp_loc / velia_study_id

        self._adata_gene = self.load_adata()
        self._adata_transcript = self.load_adata()
        
    def enable_s3(self, client:boto3.client) -> None:
        self._s3_enabled = True
        self.client = client

    def load_adatas(self) -> None:
        """
        """
        self._adata_gene = self.load_adata()
        self._adata_transcript = self.load_adata(glob_pattern='*dds_transcript*')

    def load_adata(self, glob_pattern:str='*dds_gene*') -> ad.AnnData:
        """
        """
        try:
            fh = list(Path(self.velia_study_loc / 'de_results').glob(glob_pattern))[0]
        except IndexError:
            raise FileNotFoundError
        
        if self._s3_enabled:
            response = self._client.get_object(Bucket='', Key=fh)
            response_data = response['Body'].read()
            with tempfile.TemporaryDirectory() as tempdir:
                tempfile = f'{tempdir}/gene.h5ad'
                with open(tempfile, 'wb') as f_in:
                    f_in.write(response_data)
                adata = ad.read_h5ad(tempfile, backed='r')
        else:
            adata = ad.read_h5ad(fh, backed='r')

        return adata
    
    @property
    def samples(self) -> List[str]:
        """
        """
        if (self._adata_gene.obs.index.values != self._adata_transcript.obs.index.values).all():
            raise ValueError('Study gene/transcript adatas have different samples or misordered samples.')
        
        return self._adata_gene.obs.index.to_list()

    @property 
    def samples_metadata(self) -> pd.DataFrame:
        """ 
        """
        if (self._adata_gene.obs.index.values != self._adata_transcript.obs.index.values).all():
            raise ValueError('Study gene/transcript adatas have different samples or misordered samples.')
        
        return self._adata_gene.obs.copy()

    def prepare_differentialexpression(
                self,
                de_columns:List[str]=[
                                    'baseMean','log2FoldChange','lfcSE','stat','pvalue',
                                    'padj','-log10_padj','control_mean','case_mean'],
                measurement_type:str='gene',
                ) -> Dict[str,Tuple[pd.DataFrame,List[str],List[str]]]:
        """
        """
        if measurement_type == 'transcript':
            adata = self._adata_transcript
        else:
            adata = self._adata_gene

        contrast_dict = {}
        contrasts = [c for c in adata.uns['contrasts'].keys()]
        for c in contrasts:
            use_columns = de_columns.copy()
            if 'control_mean' in de_columns and 'case_mean' in de_columns:
                try:
                    control_loc, case_loc = np.where([c.endswith('meannormedcounts') for c in adata.uns['stat_results'][c]])[0]
                    use_columns[use_columns.index('control_mean')] = adata.uns['stat_results'][c].columns[control_loc]
                    use_columns[use_columns.index('case_mean')] = adata.uns['stat_results'][c].columns[case_loc]
                except:
                    raise Exception('DE columns not formatted properly.')
            else:
                use_columns = [c for c in use_columns if c not in ('control_mean','case_mean',)]
                de_columns = [c for c in use_columns if c not in ('control_mean','case_mean',)]

            de_df = adata.uns['stat_results'][c][use_columns].copy()

            de_df.columns = [c.lower().strip('-') for c in de_columns]
            de_df['sample_condition_key'] = adata.uns['contrasts'][c][0]
            de_df['left_condition_level'] = adata.uns['contrasts'][c][1]
            de_df['right_condition_level'] = adata.uns['contrasts'][c][2]
            de_df['left_condition_display'] = c.upper().split('_VS_')[0]
            de_df['right_condition_display'] = c.upper().split('_VS_')[1]

            de_df.loc[de_df['padj'].isna(),'padj'] = 1.0
            de_df.loc[de_df['pvalue'].isna(),'pvalue'] = 1.0
            de_df.loc[de_df['log2foldchange'].isna(),'log2foldchange'] = 0.0

            left_samples = adata.obs[adata.obs[adata.uns['contrasts'][c][0]] == adata.uns['contrasts'][c][1]].index.tolist()
            right_samples = adata.obs[adata.obs[adata.uns['contrasts'][c][0]] == adata.uns['contrasts'][c][2]].index.tolist()

            contrast_dict[c.upper()] = (de_df, left_samples, right_samples)

        return contrast_dict

    def prepare_measurements(
                self, 
                measurements:List[str]=['counts','normed_counts','raw_tpm','normed_counts_transform'],
                measurement_type:str='gene',
                ) -> Tuple[np.ndarray,np.ndarray,np.ndarray,List[str]]:
        """
        """
        if measurement_type == 'transcript':
            adata = self._adata_transcript
        else:
            adata = self._adata_gene

        sequenceregions = adata.var.index.values.copy()
        samples = adata.obs.index.values.copy()
        array = np.dstack([adata.layers[l] for l in measurements])

        return np.repeat(samples, sequenceregions.shape[0]), \
                np.repeat(np.reshape(sequenceregions, (1,-1)), samples.shape[0], axis=0).flatten(), \
                np.reshape(array, (-1,len(measurements))), \
                measurements

class MetaDataFetcher:    
    """
    """
    _search_sra_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?term={srp_id}&db=sra'
    _link_sra_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?id={srp_id}&db=bioproject&dbfrom=sra'
    _link_geo_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?term={gse_id}&db=bioproject'
    _summary_bioproject_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?id={prj_id}&db=bioproject'
    _bioproject_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?id={bio_id}&db=bioproject&retmode=xml'
    _sra_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?id={srx_ids}&db=sra&rettype=runinfo&retmode=text'
    _ena_url = 'https://www.ebi.ac.uk/ena/portal/api/filereport?result=read_run&accession={bio_id}&fields={extra_params}'

    _ena_fields = [
        "run_accession",
        "experiment_accession",
        "sample_accession",
        "secondary_sample_accession",
        "study_accession",
        "secondary_study_accession",
        "submission_accession",
        "run_alias",
        "experiment_alias",
        "sample_alias",
        "library_layout",
        "library_selection",
        "library_source",
        "library_strategy",
        "library_name",
        "instrument_model",
        "instrument_platform",
        "base_count",
        "read_count",
        "tax_id",
        "scientific_name",
        "sample_title",
        "sample_description",
        ]

    _geo_id = None
    _srp_id = None
    _project_id = None
    _project_title = None
    _project_summary = None
    _sra_df = None
    _pmids = []

    def __init__(self, study_id:str, srx_ids:List[str]):
        self._study_id = study_id
        self._srx_ids = srx_ids
        self.resolve_all_ids()
    
    @property
    def samples_metadata(self) -> pd.DataFrame:
        """ 
        """
        columns_sum = ['spots', 'bases', 'spots_with_mates', 'size_MB', 'base_count', 'read_count']
        columns_mean = ['avgLength', 'InsertSize', 'InsertDev']

        df = self._sra_df.copy()

        for c in df.columns[~df.columns.isin(columns_sum+columns_mean)]:
            df[c] = df[c].astype(str)

        for c in columns_sum:
            if c not in df.columns:
                continue
            df[c] = df.groupby('Experiment')[c].transform('sum')
        for c in columns_mean:
            if c not in df.columns:
                continue
            df[c] = df.groupby('Experiment')[c].transform('mean')

        df.drop_duplicates('Experiment', inplace=True)

        return df


    def fetch_url(self, url:str, attempt_n:int=0, max_attempts:int=5) -> None:
        """ 
        """
        try:
            response = urllib.request.urlopen(url)
        except urllib.error.HTTPError as e:
            # logging.info('')
            if 'Retry-After' in e.headers:
                time.sleep(int(e.headers['Retry-After']))
            else:
                time.sleep(15)
            if attempt_n > max_attempts:
                raise Exception(f'Max attempts on url: {url}')
            attempt = attempt_n+1
            return self.fetch_url(url, attempt_n=attempt)
        return response

    def fetch_bioproject_info(self) -> http.client.HTTPResponse:
        """ 
        """
        try:
            response = self.fetch_url(self._bioproject_url.format(bio_id=self._project_id)).read()
            tree = ET.fromstring(response)
            if tree.find('.//error'):
                raise ValueError('Cannot find project_id specified.')
        except ValueError as e:
            raise e
        
        try:
            _geo_id = tree.find('.//CenterID').text
            if self._geo_id and self._geo_id != _geo_id:
                raise ValueError('Provided geo_id not the same as geo_id associated with project.')
            elif not self._geo_id:
                self._geo_id = _geo_id                
        except AttributeError as e:
            pass
        except ValueError as e:
            pass

        try:
            self._project_title = tree.find('.//ProjectDescr/Name').text
            self._project_summary = tree.find('.//ProjectDescr/Description').text
        except AttributeError as e:
            pass
    
        try:
            for pmid in tree.findall('.//Publication/Reference'):
                self._pmids.append(pmid.text)
        except AttributeError as e:
            pass

    def link_project_id(self, is_sra:bool=True) -> str:
        """ 
        """
        if is_sra:
            try:
                # Need to get an SRR id from the SRP.
                response = self.fetch_url(self._search_sra_url.format(srp_id=self._study_id)).read()
                tree = ET.fromstring(response)
                srr_id = tree.find('.//IdList/Id').text
                response = self.fetch_url(self._link_sra_url.format(srp_id=srr_id)).read()
                tree = ET.fromstring(response)
                prj_id = tree.find('.//LinkSetDb[LinkName="sra_bioproject"]/Link/Id').text
                if not prj_id:
                    raise ValueError('Unable to find project_id in bioprojects.')
            except AttributeError as e:
                raise e
            except ValueError as e:
                raise e

        else:
            try:
                response = self.fetch_url(self._link_geo_url.format(gse_id=self._study_id)).read()
                tree = ET.fromstring(response)
                prj_id = tree.find('.//IdList/Id').text
                if not prj_id:
                    raise ValueError('Unable to find project_id in bioprojects.')
            except AttributeError as e:
                raise e
            except ValueError as e:
                raise e
            
        try:
            response = self.fetch_url(self._summary_bioproject_url.format(prj_id=prj_id)).read()
            tree = ET.fromstring(response)
            bioproject_id = tree.find('.//Project_Acc').text
            if not bioproject_id:
                raise ValueError('Unable to bioproject summary from bioproject id.')
        except AttributeError as e:
            raise e
        except ValueError as e:
            raise e

        return bioproject_id

    def fetch_srx_info(self) -> None:
        """ 
        """
        try:
            response = self.fetch_url(self._sra_url.format(srx_ids=','.join(self._srx_ids))).read()
            sra_df = pd.read_csv(StringIO(response.decode()), quotechar='"', delimiter=',')
            if sra_df['Experiment'].unique().shape[0] != len(self._srx_ids):
                raise ValueError('Number of srx_ids retrieved does not equal number initialized. Probably misformatted srx_id.')
        except ValueError as e:
            raise e
        if self._study_id.startswith('GS'):
            self._srp_id = sra_df['SRAStudy'].unique()[0]

        try:
            response = self.fetch_url(self._ena_url.format(bio_id=self._project_id,extra_params=','.join(self._ena_fields))).read()
            ena_df = pd.read_csv(StringIO(response.decode()), delimiter='\t')
        except urllib.error.HTTPError as e:
            raise e

        self._sra_df = sra_df.merge(ena_df, left_on='Run', right_on='run_accession', how='left')
        
    def resolve_all_ids(self) -> None:
        """
        """
        try:
            if self._study_id.startswith('GS'):
                self._geo_id = self._study_id
                self._project_id = self.link_project_id(is_sra=False)
            elif any(map(lambda x: self._study_id.startswith(x), ('ER','SR','DR',))):
                self._srp_id = self._study_id
                self._project_id = self.link_project_id(is_sra=True)
            else:
                raise ValueError('Cannot define db from study_id.')
        except ValueError as e:
            raise e
        
        self.fetch_bioproject_info()
        self.fetch_srx_info()


if __name__ == '__main__':
    pass



